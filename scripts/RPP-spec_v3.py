#!/usr/bin/env python
from xspec import *

import matplotlib.pyplot as plt
import argparse
import os, sys
import time
import numpy as np
import re
import traceback
import copy

# -------------------- Helpers & parsing --------------------

def list_of_floats(a):
    """Parse a comma-separated string into a list of floats."""
    return list(map(float, a.split(",")))

def list_of_strings(a):
    """Parse a comma-separated string into a list of strings."""
    return a.split(",")

def insert_gaussians_in_expr(expr: str, n_gauss: int) -> str:
    """
    Build the true XSPEC expression by inserting n 'gaussian' terms.
    Insert before the last ')' if present, otherwise append at the end.
    """
    if n_gauss <= 0:
        return expr
    add = " + " + " + ".join(["gaussian"] * n_gauss)
    idx = expr.rfind(")")
    if idx != -1:
        return expr[:idx] + add + expr[idx:]
    else:
        return expr + add

def pretty_expr_with_gauss(expr: str, n_gauss: int) -> str:
    """
    Build a compact, human-readable model string for figure titles.
    Shows '+ Ngaussian' instead of repeating 'gaussian' N times.
    Mirrors the same insertion position used by insert_gaussians_in_expr.
    """
    if n_gauss <= 0:
        return expr
    add = f" + {n_gauss}gaussian"
    idx = expr.rfind(")")
    if idx != -1:
        return expr[:idx] + add + expr[idx:]
    else:
        return expr + add

def flatten(list_of_lists):
    """Flatten one level of nesting."""
    out = []
    for L in list_of_lists:
        out.extend(L)
    return out

def component_attr_names(model):
    """
    Return component attribute names in model order so that repeated components
    are accessible uniquely: e.g., 'powerlaw', 'gaussian', 'gaussian_2', ...
    """
    seen = {}
    attrs = []
    for cname in model.componentNames:
        base = cname
        seen[base] = seen.get(base, 0) + 1
        attr = base if seen[base] == 1 else f"{base}_{seen[base]}"
        attrs.append(attr)
    return attrs

def gaussian_attr_names(model):
    """Return gaussian attribute names in model order (e.g., 'gaussian', 'gaussian_2', ...)."""
    return [attr for attr in component_attr_names(model) if attr.lower().startswith("gaussian")]

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="""Do RPP catalog spectral analysis (extended with Gaussians via CLI).""",
)

parser.add_argument("model_str", type=str,
                    help="Base model string, e.g. 'tbabs*(bbody+powerlaw)'")
parser.add_argument("model_startvals", type=list_of_floats,
                    help="Comma-separated starting values for the base model (no spaces).")
parser.add_argument("--add-gauss", dest="add_gauss", action="append", default=[],
                    type=list_of_floats,
                    help="Add one Gaussian (repeatable). Format: LineE_keV,Sigma_keV,Norm")
parser.add_argument("--gauss-e-window", type=float, default=0.05,
                    help="Symmetric window ±ΔE (keV) to constrain each Gaussian LineE around its initial value.")
parser.add_argument("--gauss-sigma-bounds", type=list_of_floats, default=[0.001, 0.1],
                    help="Constrain Gaussian Sigma for all lines: 'min,max' in keV (hard+soft bounds). Default: 0.001,0.1 keV.")
parser.add_argument("--errors", action="store_true", default=False,
                    help="Run XSPEC Fit.error on Gaussian parameters as well.")
parser.add_argument("--err-max", type=int, default=100,
                    help="Maximum steps for Fit.error (used with --errors).")
parser.add_argument("--loadfile", default="merged_cutmpu7_loadRPP.py",
                    help="Python load script produced by nicerl3-spect.")
parser.add_argument("--outname", default="rpp-scorp",
                    help="Basename for output files.")
parser.add_argument("--srcname", default="",
                    help="Source name shown on the plot.")
parser.add_argument("--freeze", type=list_of_strings, default="nxb(7)",
                    help="Comma-separated parameters to freeze (no spaces), e.g. 'nxb(7),m1(12)'.")
parser.add_argument("--thaw", type=list_of_strings, default="sky(12)",
                    help="Comma-separated parameters to thaw (no spaces).")
parser.add_argument("--warn", action="store_true", default=False,
                    help="Ask before overwriting .xcm (otherwise overwrite silently).")
parser.add_argument("--timestamp", action="store_true", default=False,
                    help="Append date/time to output filenames.")

parser.add_argument("--nxb-set", action="append", default=[],
                    help="Set NXB parameter values (XSPEC syntax). Multiple assignments may be separated by ';'.")
parser.add_argument("--sky-set", action="append", default=[],
                    help="Set SKY parameter values (XSPEC syntax). Multiple assignments may be separated by ';'.")

args = parser.parse_args()

# -------------------- XSPEC globals --------------------

Plot.splashPage = False
Xset.chatter = 10
Xset.abund = "wilm"

# Load data
exec(open(args.loadfile).read())
spec1.ignore("bad")
spec1.ignore("0.0-0.25")

# Build model (true expression) and a compact string for titles
base_expr = args.model_str
n_gauss = len(args.add_gauss)
full_expr = insert_gaussians_in_expr(base_expr, n_gauss)
pretty_expr = pretty_expr_with_gauss(base_expr, n_gauss)

# Starting values = base + each Gaussian (3 params each)
startvals = args.model_startvals[:] + flatten(args.add_gauss)

m1 = Model(full_expr, modName="m1")
m1.setPars(startvals)

# ---- Remap Gaussians to CLI values + enforce hard/soft bounds on LineE and Sigma
if n_gauss > 0:
    gauss_attrs = gaussian_attr_names(m1)
    if len(gauss_attrs) != n_gauss:
        print(f"[WARN] Found {len(gauss_attrs)} Gaussian components in model "
              f"but received {n_gauss} --add-gauss entries.")

    # Parse and sanitize Sigma bounds
    if len(args.gauss_sigma_bounds) != 2:
        raise ValueError("--gauss-sigma-bounds requires 'min,max' (two floats in keV)")
    smin_cli, smax_cli = sorted(args.gauss_sigma_bounds)
    smin_cli = max(0.0, smin_cli)
    sigma_bounds = (smin_cli, smax_cli)

    for i, attr in enumerate(gauss_attrs):
        if i >= len(args.add_gauss):
            break
        comp = getattr(m1, attr)
        E0, S0, N0 = args.add_gauss[i]

        # LineE: clamp within ±window around E0 (both hard and soft limits)
        if args.gauss_e_window and args.gauss_e_window > 0:
            vE = list(comp.LineE.values)  # [val, delta, min, bot, top, max]
            vE[0] = E0
            new_min = max(E0 - args.gauss_e_window, vE[2])
            new_max = min(E0 + args.gauss_e_window, vE[5])
            if new_min > new_max:
                new_min, new_max = new_max, new_min
            vE[2] = vE[3] = new_min
            vE[4] = vE[5] = new_max
            vE[1] = max(args.gauss_e_window / 5.0, 1e-4) if (vE[1] <= 0 or vE[1] > args.gauss_e_window / 5.0) else vE[1]
            comp.LineE.values = ",".join(str(x) for x in vE)
        else:
            comp.LineE.values = f"{E0},{max(E0,1e-4)/10.0}"

        # Sigma: enforce global bounds (both hard and soft) and clip start
        vS = list(comp.Sigma.values)
        smin, smax = sigma_bounds
        smin_eff = max(smin, vS[2])
        smax_eff = min(smax, vS[5])
        if smin_eff > smax_eff:
            smin_eff, smax_eff = smax_eff, smin_eff
        vS[2] = vS[3] = smin_eff
        vS[4] = vS[5] = smax_eff
        S0 = min(max(S0, smin_eff), smax_eff)
        vS[0] = S0
        if vS[1] <= 0:
            vS[1] = max(S0, 1e-3) / 10.0
        comp.Sigma.values = ",".join(str(x) for x in vS)

        # Norm: set start value and a reasonable positive step
        vN = list(comp.norm.values)
        vN[0] = N0
        if vN[1] <= 0:
            vN[1] = max(abs(N0), 1e-8) / 10.0
        comp.norm.values = ",".join(str(x) for x in vN)

AllModels.show()

# Separate handles for sky and NXB
sky = AllModels(1, modName="sky")
nxb = AllModels(1, modName="nxb")


timestamp = time.strftime("%Y%m%d-%H%M%S")

# -------------------- Fit --------------------

if args.timestamp:
    logfile = f"{args.outname}-{timestamp}.log"
    xcmfile = f"{args.outname}-{timestamp}.xcm"
    ymlfile = f"{args.outname}-{timestamp}.yml"
    pngfile = f"{args.outname}-{timestamp}.png"
else:
    logfile = f"{args.outname}.log"
    xcmfile = f"{args.outname}.xcm"
    ymlfile = f"{args.outname}.yml"
    pngfile = f"{args.outname}.png"

frozen_list = []

# --- Background presets (values + keep them thawed) -----------------------
import re as _re

def _expand_assignments(list_or_none):
    """Split each provided string on ';' and strip."""
    out = []
    for s in (list_or_none or []):
        parts = [p.strip() for p in s.split(";") if p.strip()]
        out.extend(parts)
    return out

def _apply_assignments(assign_list, ctx_name, modified_params, frozen_log):
    """
    Each assignment is like:
      'nxb(4)=8.7'  or  'nxb(4)=8.7,0.05,0,0,1e20,1e24'  or  'sky(6)=0'
    LHS is evaluated with eval() (e.g. nxb(4)); RHS is passed to Parameter.values.
    The parameter is made thawed (frozen=False).
    """
    for expr in _expand_assignments(assign_list):
        try:
            lhs, rhs = _re.split(r"\s*=\s*|\s*:\s*", expr, maxsplit=1)
            lhs = lhs.strip()
            rhs = rhs.strip()
            p = eval(lhs)              # e.g. nxb(4) or sky(6)
            p.values = rhs             # can be 'val' or 'val,delta,min,bot,top,max'
            p.frozen = False           # ensure this is free
            modified_params.append(p)  # remember for post-FREEZE re-thaw
            tmp = " ".join(map(str, p.values))
            frozen_log.append(f"\n[PRESET {ctx_name}] {lhs} set -> (val,delta,min,bot,top,max): {tmp}  [thawed]")
        except Exception as e:
            frozen_log.append(f"\n[WARN PRESET {ctx_name}] Could not apply '{expr}': {e}")

# collect modified background params to re-thaw them later if needed
_modified_bg_params = []
_apply_assignments(args.nxb_set, "NXB", _modified_bg_params, frozen_list)
_apply_assignments(args.sky_set, "SKY", _modified_bg_params, frozen_list)
# -------------------------------------------------------------------------






# Pre-fit THAW
frozen_list.append("\n\nTHAWED on command line by user pre-fit:\n")
for pname in (args.thaw if isinstance(args.thaw, list) else [args.thaw]):
    if not pname:
        continue
    p = eval(pname)
    p.frozen = False
    tmp = " ".join(list(map(str, p.values)))
    frozen_list.append(f"\n{pname} thawed with (val,delta,min,bot,top,max): {tmp}")

# Pre-fit FREEZE
frozen_list.append("\n\nFROZEN on command line by user pre-fit:\n")
for pname in (args.freeze if isinstance(args.freeze, list) else [args.freeze]):
    if not pname:
        continue
    p = eval(pname)
    p.frozen = True
    tmp = " ".join(list(map(str, p.values)))
    frozen_list.append(f"\n{pname} frozen with (val,delta,min,bot,top,max): {tmp}")


# Ensure background presets stay thawed even if user froze them in --freeze
for p in _modified_bg_params:
    try:
        p.frozen = False
    except Exception:
        pass




Fit.statMethod = "pgstat"
Fit.nIterations = 200
Fit.query = "no"
Fit.perform()
Fit.show()

# Conservative auto-handling of problematic additive norms only
frozen_list.append(
    "\n\nAUTO handling: set-to-zero or freeze ONLY additive 'norm' with sigma=-1 or sigma>>value:\n"
)

def is_add_norm(model_obj, par_index):
    """
    Return True if the parameter is 'norm' of a typical additive component
    (gaussian, powerlaw, bbody, bbodyrad, bknpower, apec).
    """
    p = model_obj(par_index)
    try:
        comp = p.ownerComponent
        pname = p.name.lower()
        cname = comp.name.lower()
    except Exception:
        return False
    if pname != "norm":
        return False
    return any(tag in cname for tag in ["gaussian", "powerlaw", "bbody", "bbodyrad", "bknpower", "apec"])

for m in [m1, sky, nxb]:
    touched = False
    for ii in range(m.startParIndex, m.nParameters + 1):
        p = m(ii)
        if not is_add_norm(m, ii):
            continue  # do not touch LineE/Sigma/etc.
        if (p.sigma == -1) or (p.values[0] != 0 and p.sigma > 10.0 * abs(p.values[0])):
            if p.values[3] <= 0 and p.values[5] >= 0:
                p.values = "0.0,0"
            else:
                p.frozen = True
            tmp = " ".join(list(map(str, p.values)))
            frozen_list.append(f"\n{m.name}:{p.ownerComponent.name}.{p.name} handled -> (val,delta,min,bot,top,max): {tmp}")
            touched = True
    if touched:
        Fit.perform()

# Optional: XSPEC errors for Gaussian parameters (using indexed attrs)
if args.errors:
    gauss_attrs = gaussian_attr_names(m1)
    gauss_par_ids = []
    for attr in gauss_attrs:
        c = getattr(m1, attr)
        for pn in ("LineE", "Sigma", "norm"):
            p = getattr(c, pn)
            if not p.frozen:
                gauss_par_ids.append(str(p.index))
    if gauss_par_ids:
        Fit.error("maximum {} {}".format(args.err_max, " ".join(gauss_par_ids)))

# Log + flux
Xset.openLog(logfile)
AllData.show()
AllModels.show()
Fit.show()
AllModels.calcFlux("0.4 2.0")
flux_loint = AllData(1).flux[0]
AllModels.calcFlux("2.0 10.0")
flux_hiint = AllData(1).flux[0]
Xset.closeLog()
print("Log file closed.")

with open(logfile, "a+") as fd:
    for f in frozen_list:
        fd.write(f)

if os.path.exists(xcmfile) and not args.warn:
    os.remove(xcmfile)
Xset.save(xcmfile)

# -------------------- Save YAML --------------------
# Use component attribute names so that repeated components are handled correctly.

bestfit_params = []

with open(ymlfile, "w") as fd:
    bbody_count = 0
    gauss_count = 0

    comp_attrs = component_attr_names(m1)
    for attr in comp_attrs:
        c = getattr(m1, attr)
        base = attr.split("_")[0].lower()

        if base == "bbody":
            bbody_count += 1
        if base == "gaussian":
            gauss_count += 1

        for pname in c.parameterNames:
            p = getattr(c, pname)
            pval = p.values[0]
            psigma = p.sigma
            plo = max(p.values[0] - p.sigma, p.values[2])
            phi = min(p.values[0] + p.sigma, p.values[5])

            bestfit_params.append(pval)
            
            if base == "tbabs":
                fd.write("nH: %g\n" % pval)
                fd.write("nH_sigma: %g\n" % psigma)
                fd.write("nH_high: %g\n" % phi)
                fd.write("nH_low: %g\n" % plo)

            elif base == "bbody":
                # First bbody uses suffix '', others use index
                if bbody_count < 2:
                    fd.write("%s_thermal: %g\n" % (pname, pval))
                    fd.write("%s_thermal_sigma: %g\n" % (pname, psigma))
                    fd.write("%s_thermal_high: %g\n" % (pname, phi))
                    fd.write("%s_thermal_low: %g\n" % (pname, plo))
                else:
                    fd.write("%s_thermal%d: %g\n" % (pname, bbody_count, pval))
                    fd.write("%s_thermal%d_sigma: %g\n" % (pname, bbody_count, psigma))
                    fd.write("%s_thermal%d_high: %g\n" % (pname, bbody_count, phi))
                    fd.write("%s_thermal%d_low: %g\n" % (pname, bbody_count, plo))

            elif base == "powerlaw":
                pn = "gamma" if "PhoIndex" in pname else pname
                fd.write("%s_powerlaw: %g\n" % (pn, pval))
                fd.write("%s_powerlaw_sigma: %g\n" % (pn, psigma))
                fd.write("%s_powerlaw_high: %g\n" % (pn, phi))
                fd.write("%s_powerlaw_low: %g\n" % (pn, plo))

            elif base == "gaussian":
                fd.write("gauss%d_%s: %g\n" % (gauss_count, pname, pval))
                fd.write("gauss%d_%s_sigma: %g\n" % (gauss_count, pname, psigma))
                fd.write("gauss%d_%s_high: %g\n" % (gauss_count, pname, phi))
                fd.write("gauss%d_%s_low: %g\n" % (gauss_count, pname, plo))

    fd.write("flux_0.4-2: %g\n" % flux_loint)
    fd.write("flux_2-10: %g\n" % flux_hiint)

# -------------------- Plot --------------------

# Pull arrays from XSPEC (no on-screen XSPEC plotting)
Plot.device = "/null"
Plot.xAxis = "keV"
Plot.add = True
Plot.xLog = True
Plot.yLog = True
Plot.setRebin()
Plot.addCommand("res y 0.01 20")
Plot("ldata")
en = Plot.x()
rates = Plot.y()
rates_err = Plot.yErr()
folded = Plot.model()

# Ratio arrays
Plot.yLog = False
Plot("ratio")
en_ratio = Plot.x()
ratio = Plot.y()
ratio_err = Plot.yErr()

# Bin widths (for counts/s/keV in component plots)
en = np.array(en)
en_dif = en[1:] - en[:-1]
en_dif = np.concatenate((en_dif, np.array([en_dif[-1]])))

fig, (ax1, ax2) = plt.subplots(
    2, 1, sharex=True, gridspec_kw={"hspace": 0, "height_ratios": [3, 1]}
)
# ---- Use compact model string in the figure title
fig.suptitle(f"{args.srcname}   Model: {pretty_expr}")

ax1.errorbar(en, rates, yerr=rates_err, color="black", linestyle="")
ax1.step(en, rates, label="Data", color="black")
ax1.plot(en, folded, label="Full Model", color="orange")
ax1.plot(en, np.divide(m1.folded(1), en_dif), label="Source", color="red")
ax1.plot(en, np.divide(sky.folded(1), en_dif), label="Sky", color="blue")
ax1.plot(en, np.divide(nxb.folded(1), en_dif), label="NXB", color="green")
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_ylim(bottom=1.0e-6)
ax1.tick_params(axis="x", which="both", bottom=True, top=True, direction="in", labelbottom=False)

# Plot each Gaussian by zeroing other additive norms (use indexed attrs)
if n_gauss > 0:
    comp_attrs = component_attr_names(m1)
    comps = []
    for attr in comp_attrs:
        comp_obj = getattr(m1, attr)
        pnorm = getattr(comp_obj, "norm") if "norm" in comp_obj.parameterNames else None
        comps.append((attr.lower(), comp_obj, pnorm))

    additive_norm_params = [(nm, p) for (nm, comp, p) in comps if p is not None]
    gauss_attrs = gaussian_attr_names(m1)
    gauss_norms = [getattr(getattr(m1, attr), "norm") for attr in gauss_attrs]

    for gidx, gnorm in enumerate(gauss_norms, start=1):
        saved_vals = [p.values[0] for (_, p) in additive_norm_params]
        for _, p in additive_norm_params:
            p.values = "0.0,0"
        try:
            gi = [p for (_, p) in additive_norm_params].index(gnorm)
            gsv = saved_vals[gi]
        except ValueError:
            gsv = gnorm.values[0]
        gnorm.values = f"{gsv},0"
        line_only = np.divide(m1.folded(1), en_dif)
        ax1.plot(en, line_only, linestyle="--", linewidth=1.0, label=f"Gauss {gidx}")
        for (nm, p), sv in zip(additive_norm_params, saved_vals):
            p.values = f"{sv},0"

ax1.legend(ncol=5, mode="expand")
ax1.set_ylabel("counts / s / keV")
ax1.set_title(
    "Chi2: %.2f  DOF: %d  Chi2/DOF: %.2f  S[0.4-2]: %.2g  S[2-10]: %.2g"
    % (Fit.testStatistic, Fit.dof, Fit.testStatistic / Fit.dof, flux_loint, flux_hiint),
    fontsize="medium",
)

ax2.errorbar(en_ratio, ratio, yerr=ratio_err, color="black")
ax2.axhline(1, color="blue", linewidth=0.5)
ax2.set_xscale("log")
ax2.set_ylim(0.75, 1.25)
ax2.set_xlabel("Energy (keV)")
ax2.set_ylabel("Data/Full Model")
ax2.tick_params(axis="x", which="both", top=True, direction="in", labelbottom=True)

plt.savefig(pngfile)


# cflux calculations -------------------------------------------------------

# Args: the model string WITH cflux; the best-fit model parameters from the previously done fit of the model WITHOUT cflux; cflux starting parameters
def cflux(cflux_model_str,model_pars,cflux_start_pars,loadfile=args.loadfile,models=AllModels):
    print('\n')
    print('**************************************************')
    print('****************CFLUX CALCULATION*****************')
    print('**************************************************')

    # Figure out where in the starting model parameters array we need to insert the cflux starting params
    comps = re.split('[^a-zA-Z]', cflux_model_str)
    count_pars = 0 
    for comp in comps:
        if 'tbabs' in comp.lower():
            count_pars += 1
        elif 'body' in comp.lower() or 'pow' in comp.lower():
            count_pars += 2
        elif 'gauss' in comp.lower():
            count_pars += 3
        elif 'cflux' in comp.lower():
            start_pars = model_pars[0:count_pars] + cflux_start_pars + model_pars[count_pars:]
            break

    print('Model str: ',cflux_model_str)
    print('Model start params: ',start_pars)
    print('**************************************************')
    
    # Initialize model - this assumes the loadfile has been read in the caller AFTER the plotting block. The background file is removed in the plotting block to get the src-only curve, so have to read it again before fitting the cflux models.

    # Can't do AllModels.clear() here b/c we need the sky and nxb models
    try:
        models -= "m2" #this throws an error on the 1st cflux() invocation b/c m2 doesn't yet exist; that is ok.
    except:
        #traceback.print_exc()
        pass
    
    m2 = Model(cflux_model_str, modName="m2")
    m2.setPars(start_pars)

    # Freeze Emin,Emax, and all non-cflux params (or just norms).
    # Save the flux param index - we'll use it later to get the best-fit values.
    ii_flux = -1
    for ii in range(m2.startParIndex,m2.nParameters+1):
        #if 'lg10flux' not in (m2(ii).name).lower():
        if any(elem in (m2(ii).name).lower() for elem in ['emin','emax','norm']):
            m2(ii).frozen = True
        if 'lg10flux' in (m2(ii).name).lower():
            ii_flux = ii
            #print('ii_flux: ',ii_flux)

    AllModels.show()
    
    # Do the fit
    Fit.statMethod = "pgstat"
    Fit.nIterations = 200
    Fit.query = "no"
    try:
        Fit.perform()
    except Exception as e:
        print('Fit failed.')
        #print(e)
        traceback.print_exc()
        return [-1,-1]
        
    Fit.show()
    return [m2(ii_flux).values[0], m2(ii_flux).sigma]
    

print("********************* CFLUX STARTING *****************************")

# Remove the non-cflux model, we have saved its best-fit params and don't need it any more
AllModels -= "m1"
# It is not necessary to re-load the initial file defining the data here b/c unlike the OnOff case, there is no background file removal in the plotting block. All data files are still in place.

# Freeze all params in the sky and nxb models. Most cases run fine if these aren't fully frozen, but I had one case go into an infinite loop.
for ii in range(nxb.startParIndex,nxb.nParameters+1):
    nxb(ii).frozen = True
for ii in range(sky.startParIndex,sky.nParameters+1):
    sky(ii).frozen = True

cflux_start_params_loint = [0.4, 2.0, np.log10(flux_loint)]
cflux_start_params_hiint = [2.0, 10.0, np.log10(flux_hiint)]
fd = open(ymlfile,'a')

# Absorbed for the full model
cflux_model_str = 'cflux*'+full_expr
log10flux_abs, log10flux_abs_sigma = cflux(cflux_model_str,bestfit_params,cflux_start_params_loint)
fd.write("log10flux_abs_0.4-2: %g\n" % log10flux_abs)
fd.write("log10flux_abs_0.4-2_sigma: %g\n" % log10flux_abs_sigma)
log10flux_abs, log10flux_abs_sigma = cflux(cflux_model_str,bestfit_params,cflux_start_params_hiint)
fd.write("log10flux_abs_2-10: %g\n" % log10flux_abs)
fd.write("log10flux_abs_2-10_sigma: %g\n" % log10flux_abs_sigma)

# Unabsorbed for the full model
cflux_model_str_unabs = 'tbabs*cflux'+(full_expr).split('tbabs')[-1]
log10flux_unabs, log10flux_unabs_sigma = cflux(cflux_model_str_unabs,bestfit_params,cflux_start_params_loint)
fd.write("log10flux_unabs_0.4-2: %g\n" % log10flux_unabs)
fd.write("log10flux_unabs_0.4-2_sigma: %g\n" % log10flux_unabs_sigma)
log10flux_unabs, log10flux_unabs_sigma = cflux(cflux_model_str_unabs,bestfit_params,cflux_start_params_hiint)
fd.write("log10flux_unabs_2-10: %g\n" % log10flux_unabs)
fd.write("log10flux_unabs_2-10_sigma: %g\n" % log10flux_unabs_sigma)

# Unabsorbed for components separately if there is more than one within braces
bbody_count = 0
gauss_count = 0

if '(' in full_expr:
    tmp1 = (full_expr).split('(')
    tmp2 = tmp1[-1].split(')')[0]
    comps = tmp2.split('+')

    for ii in range(0,len(comps)):
        newcomps = copy.deepcopy(comps)
        newcomps[ii] = 'cflux*'+comps[ii]
        newcomps = '+'.join(newcomps)

        cflux_comp_unabs = tmp1[0]+'('+newcomps+')'

        log10flux_loint, log10flux_loint_sigma = cflux(cflux_comp_unabs,bestfit_params,cflux_start_params_loint)
        log10flux_hiint, log10flux_hiint_sigma = cflux(cflux_comp_unabs,bestfit_params,cflux_start_params_hiint)

        if 'bbody' in comps[ii]:
            bbody_count += 1
            if bbody_count < 2:
                fd.write("log10flux_bbody_0.4-2: %g\n" % log10flux_loint)
                fd.write("log10flux_bbody_0.4-2_sigma: %g\n" % log10flux_loint_sigma)
                fd.write("log10flux_bbody_2-10: %g\n" % log10flux_hiint)
                fd.write("log10flux_bbody_2-10_sigma: %g\n" % log10flux_hiint_sigma)
            else:
                fd.write("log10flux_bbody%d_0.4-2: %g\n" % (bbody_count,log10flux_loint))
                fd.write("log10flux_bbody%d_0.4-2_sigma: %g\n" % (bbody_count,log10flux_loint_sigma))
                fd.write("log10flux_bbody%d_2-10: %g\n" % (bbody_count,log10flux_hiint))
                fd.write("log10flux_bbody%d_2-10_sigma: %g\n" % (bbody_count,log10flux_hiint_sigma))
        elif 'pow' in comps[ii]:
            fd.write("log10flux_pow_0.4-2: %g\n" % log10flux_loint)
            fd.write("log10flux_pow_0.4-2_sigma: %g\n" % log10flux_loint_sigma)
            fd.write("log10flux_pow_2-10: %g\n" % log10flux_hiint)
            fd.write("log10flux_pow_2-10_sigma: %g\n" % log10flux_hiint_sigma)
        elif 'gauss' in comps[ii]:
            gauss_count += 1
            fd.write("log10flux_gauss%d_0.4-2: %g\n" % (gauss_count,log10flux_loint))
            fd.write("log10flux_gauss%d_0.4-2_sigma: %g\n" % (gauss_count,log10flux_loint_sigma))
            fd.write("log10flux_gauss%d_2-10: %g\n" % (gauss_count,log10flux_hiint))
            fd.write("log10flux_gauss%d_2-10_sigma: %g\n" % (gauss_count,log10flux_hiint_sigma))
        

fd.close()

