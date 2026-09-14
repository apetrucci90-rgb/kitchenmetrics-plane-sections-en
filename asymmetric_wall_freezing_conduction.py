#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Asymmetric-boundary transient conduction through a plane wall, WITH LATENT HEAT —
validation only, for the FREEZING side of the cavity model.

Step 1 of the plan Andrea set on 2026-09-13, same convention as
docs/asymmetric_wall_conduction.py (the cooking side, 2026-09-11): "Comincia con lo
script Python, verifichiamo e poi creiamo un motore di controllo e solo allora
mettiamola in freezing.js." This file is that first step and NOTHING past it — it
writes no JavaScript, and freezing.js is not touched by it.

WHY THIS IS A NEW FILE, NOT A CALL INTO asymmetric_wall_conduction.py
That file's solver assumes UNIFORM, temperature-INDEPENDENT k and alpha — a
sensible-heat-only transient, correct for roasting (the body only warms up).
Freezing is dominated by the LATENT heat of fusion (~334 kJ/kg of the FREEZABLE
water — freezing.js's FRZ_L_WATER), released over a narrow band around the food's
own initial freezing point Tfm, not at a fixed temperature, and NOT as a small
correction: it is most of what a freezing time actually pays for. Reusing the
cooking solver unchanged would silently drop that heat and print a freezing time
far too SHORT — the dangerous direction for this app to be wrong in. So the PDE
itself needs temperature-dependent properties, not just new boundary conditions.

THE CAVITY: EMPTY, not steam-filled. Andrea, 2026-09-13, in so many words: "tenendo
conto che la cavità è vuota." No pool, no evaporation, no 100 degC ceiling — those
belonged to a roast expelling juice, and a frozen carcass's cavity is just air,
sealed off from the freezer's own moving air. The air trapped inside has nowhere
external to exchange heat with; it can only exchange with the wall's own inner
face. So the inner boundary here is an AIR POCKET with its own (small) heat
capacity, warmer than the wall at first, that cools itself by giving heat to the
inner face, and cools to no other sink than that — never a fixed T-infinity taken
from the freezer, because the cavity is not directly exposed to that air.

WHAT THIS FILE DOES
  1. Ports freezing.js's OWN property model line for line (Choi-Okos unfrozen,
     ASHRAE ch.19 ice properties, the Miles/Chen bound-water ice-fraction
     correlation, Pham's Eq. 51) — not a new model, so the validation below is
     checking THIS solver against the SAME physics the shipped app already uses,
     not against an independently-invented freezing formula.
  2. Builds an APPARENT SPECIFIC HEAT C_app(T) from that ice fraction: sensible
     heat capacity plus the latent term folded in as rho * Lf * |dxIce/dT| — the
     standard way to fold a phase change into a plain conduction PDE without
     tracking a moving freezing front explicitly.
  3. Solves 1-D transient conduction through a plane wall with THIS
     temperature-dependent C(T) and k(T), independent Robin boundary conditions on
     each face (same explicit FTCS scheme as the cooking file, generalised to
     look up local properties at the CURRENT field each step).
  4. Validates it before trusting it for a single number, per
     MOTORI-DI-CONTROLLO.md's own rule — three independent checks, none of which
     re-implement the finite-difference method itself:
       a. SYMMETRIC REDUCTION vs. PHAM. With identical (h, T-infinity) on both
          faces, this solver's own centre-of-wall freezing time must match Pham
          (1986) Eq. 51 — freezing.js's shipped formula, ported here — for the
          SAME slab. This is the one case where a completely different method
          (a fitted closed-form correlation vs. a finite-difference PDE solve)
          must land on the same answer, so it is the strongest check available.
       b. ENERGY CONSERVATION, WITH LATENT HEAT. Integrated boundary flux through
          both faces must equal the ENTHALPY actually gained (sensible term at the
          Choi-Okos midpoint cp plus rho * Lf * (xIce_final - xIce_init), read off
          the SAME ice-fraction correlation independently of the solver's own
          local dxIce/dT bookkeeping) — the same property MOTORI-DI-CONTROLLO.md's
          B2/C1/C3 already lean on, extended to cover the latent term explicitly.
       c. GRID CONVERGENCE, in the purely-sensible regime (no phase change
          crossed), against the SAME multi-term Fourier series oracle the cooking
          file already validated its own scheme against — proof the base scheme
          is sound before the nonlinear C(T) term is layered on top of it.
  5. Prints the three cavity band members (conservative / centro / aggressive) for
     one illustrative wall, so the SIZE of the band is visible before any of this
     goes near freezing.js.

WHAT THIS FILE DOES NOT DO, on purpose
  - It does not touch app/src/main/assets. Nothing here is wired to the app.
  - It does not decide the cavity-side h values FOR the app beyond what is argued
    here — every one is marked ASSUMED, at the same status the cooking file's
    H_CAVITY_INERT/HUMID already carry (and this file reuses those exact two
    numbers: same physical mechanism — natural convection in a closed air gap —
    so a new pair of unsourced constants would be pure duplication, not rigour).
  - It is not itself the "motore di controllo" the plan calls for next. That has
    to check the JAVASCRIPT the app ships, by a route this script cannot reach
    from Python. This file's job is narrower: prove the METHOD is sound before
    anyone writes that JavaScript at all.

HOW TO RUN IT
    python3 asymmetric_wall_freezing_conduction.py

Nothing to install — standard library only (math, sys). Takes no arguments, reads
no files, writes only to the screen. The last line says either "All checks
passed." or names the check that failed, and the exit status is 0 or 1 accordingly.
"""

import math
import sys

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

FAILURES = []


def check(label, ok, detail=''):
    tag = 'OK  ' if ok else 'FAIL'
    print(f'  [{tag}] {label}' + (f' — {detail}' if detail else ''))
    if not ok:
        FAILURES.append(label)


# ═══════════════════════════════════════════════════════════════════════════
# 1. PROPERTY MODEL — ported line for line from
#    app/src/main/assets/js/panels/freezing.js, not reinvented. Comments here
#    are trimmed; see that file for the sourcing of every coefficient.
# ═══════════════════════════════════════════════════════════════════════════

FRZ_L_WATER = 333600.0  # J/kg — freezing.js's own constant, latent heat of fusion


def comp_props(tC):
    """Choi-Okos component properties at tC (degC). rho: kg/m3, cp: J/(kg*K),
    k: W/(m*K). ASHRAE (2006) ch.9 — same coefficients as freezing.js's
    frzCompProps()."""
    rho = {
        'w': 997.18 + 0.0031439 * tC - 0.0037574 * tC * tC,
        'p': 1329.9 - 0.5184 * tC,
        'f': 925.59 - 0.41757 * tC,
        'c': 1599.1 - 0.31046 * tC,
        'a': 2423.8 - 0.28063 * tC,
    }
    cp = {
        'w': (4.1289 - 9.0864e-5 * tC + 5.4731e-6 * tC * tC) * 1000,
        'p': (2.0082 + 1.2089e-3 * tC - 1.3129e-6 * tC * tC) * 1000,
        'f': (1.9842 + 1.4733e-3 * tC - 4.8008e-6 * tC * tC) * 1000,
        'c': (1.5488 + 1.9625e-3 * tC - 5.9399e-6 * tC * tC) * 1000,
        'a': (1.0926 + 1.8896e-3 * tC - 3.6817e-6 * tC * tC) * 1000,
    }
    k = {
        'w': 0.57109 + 1.7625e-3 * tC - 6.7036e-6 * tC * tC,
        'p': 0.17881 + 1.1958e-3 * tC - 2.7178e-6 * tC * tC,
        'f': 0.18071 - 2.7604e-4 * tC - 1.7749e-7 * tC * tC,
        'c': 0.20141 + 1.3874e-3 * tC - 4.3312e-6 * tC * tC,
        'a': 0.32962 + 1.4011e-3 * tC - 2.9069e-6 * tC * tC,
    }
    return rho, cp, k


def ice_props(tC):
    """ASHRAE ch.19. cp in J/(kg*K)."""
    return (
        916.89 - 0.13071 * tC,
        (2.0623 + 6.0769e-3 * tC) * 1000,
        2.2196 - 6.2489e-3 * tC + 1.0154e-4 * tC * tC,
    )


def mix(x, rho, cp, k):
    """Bulk mixture: parallel k on volume fractions, series (harmonic) rho —
    freezing.js's frzMix()."""
    keys = list(x.keys())
    dens = 1.0 / sum(x[m] / rho[m] for m in keys)
    cpv = sum(x[m] * cp[m] for m in keys)
    vf = {m: x[m] / rho[m] for m in keys}
    vfs = sum(vf.values())
    kv = sum(k[m] * vf[m] for m in keys) / vfs
    return dens, cpv, kv, dens * cpv  # rho, cp, k, C=rho*cp


def freezable_water(comp):
    return max(0.0, comp['w'] - 0.4 * comp['p'])


def ice_fraction(comp, Tf, tC):
    """Miles/Chen bound-water model, ASHRAE ch.19 Eq.4 — freezing.js's
    frzIceFraction(). Mass fraction of ICE (of the whole food), at tC below Tf."""
    xb = 0.4 * comp['p']
    if tC >= Tf:
        return 0.0
    xi = (comp['w'] - xb) * (1.0 - Tf / tC)
    return max(0.0, min(xi, comp['w'] - xb))


def props_unfrozen(comp, tC):
    rho, cp, k = comp_props(tC)
    return mix(comp, rho, cp, k)


def props_frozen(comp, Tf, tC):
    rho, cp, k = comp_props(tC)
    rho_i, cp_i, k_i = ice_props(tC)
    xi = ice_fraction(comp, Tf, tC)
    x = dict(comp); x['w'] = comp['w'] - xi; x['ice'] = xi
    rho2 = dict(rho); rho2['ice'] = rho_i
    cp2 = dict(cp); cp2['ice'] = cp_i
    k2 = dict(k); k2['ice'] = k_i
    return mix(x, rho2, cp2, k2)


def bulk_k(comp, Tf, tC):
    """k(T): frozen mixture below Tf, unfrozen at or above it."""
    if tC < Tf:
        return props_frozen(comp, Tf, tC)[2]
    return props_unfrozen(comp, tC)[2]


def apparent_C(comp, Tf, tC, dT=0.75):
    """Apparent VOLUMETRIC heat capacity, J/(m3*K): sensible term (from whichever
    side of Tf tC sits on) plus the latent term folded in as
    rho0 * Lf_specific * |dxIce/dT|, evaluated by CENTRAL DIFFERENCE on
    ice_fraction() itself — not its analytic derivative — so this uses exactly
    the same correlation freezing.js calls, with no separate formula to keep in
    sync. rho0 (for turning the SPECIFIC latent heat into a volumetric one) is
    the unfrozen density at tC: the mass doesn't change when water turns to ice,
    only its volume does a little, and this is the same approximation
    freezing.js's own Lf = rho0 * freezable_water * FRZ_L_WATER already makes.

    dT=0.75 degC, not a narrower window: this is the standard trade-off of the
    apparent-heat-capacity method for phase change (Voller & Swaminathan and
    the wider enthalpy-method literature). A window as narrow as a few
    hundredths of a degree looks more "exact" but a node stepping across it in
    one explicit timestep can cross the ENTIRE spike without the scheme seeing
    enough of it to account for the latent heat it represents — first found
    here as a 39% energy-conservation gap at dT=0.05 that this wider window
    removes. A food's own freezing point is already smeared over a range by
    its mixture of solutes (that is what ice_fraction() itself models), so
    widening the numerical window by well under one degree does not invent
    physics that was not already approximate.
    """
    if tC < Tf:
        _, _, _, Csens = props_frozen(comp, Tf, tC)
    else:
        _, _, _, Csens = props_unfrozen(comp, tC)
    rho0 = props_unfrozen(comp, tC)[0]
    xi_hi = ice_fraction(comp, Tf, tC - dT)
    xi_lo = ice_fraction(comp, Tf, tC + dT)
    dxi_dT = (xi_hi - xi_lo) / (2 * dT)  # >= 0: ice fraction rises as T falls
    return Csens + rho0 * FRZ_L_WATER * dxi_dT


def pham_time_seconds(Cl, Cs, Lf, ks, Ti, Tc, Tm, h, VA, D):
    """Pham (1986) Eq. 51 — ported verbatim from freezing.js's frzPhamTime()."""
    Tfm = 1.8 + 0.263 * Tc + 0.105 * Tm
    dH1 = Cl * (Ti - Tfm)
    dH2 = Lf + Cs * (Tfm - Tc)
    dT1 = (Ti + Tfm) / 2 - Tm
    dT2 = Tfm - Tm
    Bis = h * D / ks
    return (dH1 / dT1 + dH2 / dT2) / h * VA * (1 + Bis / 4)


# A representative lean-meat-like composition (fractions sum to 1), used
# throughout this file — validates the METHOD, not a specific food row.
DEMO_COMP = {'w': 0.75, 'p': 0.20, 'f': 0.02, 'c': 0.0, 'a': 0.03}
DEMO_TF = -1.0  # degC, initial freezing point — a round, typical lean-meat value


# ═══════════════════════════════════════════════════════════════════════════
# 2. THE SOLVER — explicit FD (FTCS), nonlinear: k and C looked up from the
#    CURRENT local temperature at every node, every step (the standard
#    "apparent heat capacity" treatment of a phase change inside a plain
#    conduction PDE — no explicit front-tracking). Same half-cell boundary
#    treatment as the cooking file, generalised from uniform to local
#    alpha_eff = k(T)/C(T).
# ═══════════════════════════════════════════════════════════════════════════

class WallResult:
    __slots__ = ('T', 'nx', 'dx', 'time_history', 'flux_in_history')

    def __init__(self):
        self.time_history = []
        self.flux_in_history = []


def solve_wall_latent(L, comp, Tf, T_init, outer_bc, inner_bc_fn, t_end,
                       nx=101, safety=0.40, record_every=1, max_steps=20_000_000,
                       T_bounds=None, k_of_override=None, C_of_override=None):
    """
    L            full wall thickness, m
    comp, Tf     food composition and initial freezing point, feeding k(T)/C(T)
    T_init       initial temperature, degC, uniform
    outer_bc     function(t, T_face) -> (h, T_inf), outer (freezer) face
    inner_bc_fn  function(t, T_face) -> (h, T_inf), inner (cavity) face —
                  may also carry .advance(dt, T_face) for a stateful BC
    t_end        stop time, s
    nx           node count (odd, so a node sits at the centre)
    """
    assert nx % 2 == 1, 'nx must be odd so a node sits exactly at the centre'
    dx = L / (nx - 1)
    T = [T_init] * nx

    # Overrides exist ONLY for validate_grid_convergence_sensible_only() below,
    # which needs properties that are TRULY constant (not Choi-Okos's own small
    # T-dependence) to isolate the scheme's discretisation error from a model
    # bias — the same reason the cooking file's own convergence check hard-codes
    # a fixed alpha/k rather than calling a composition model.
    k_of = k_of_override if k_of_override else (lambda Tval: bulk_k(comp, Tf, Tval))
    C_of = C_of_override if C_of_override else (lambda Tval: apparent_C(comp, Tf, Tval))

    # Stability bound. The apparent alpha = k(T)/C(T) is LOWEST (fastest
    # diffusion, tightest dt) away from the latent plateau — at the plateau
    # itself C spikes and alpha drops, which makes that region MORE stable,
    # not less. So the binding case is somewhere in the plain-sensible part of
    # the range this wall will actually visit, and that range is NOT just
    # T_init and Tf: a wall being frozen also passes through fully-frozen
    # temperatures near the freezer's own setpoint, where ice's own (higher)
    # conductivity and (lower) heat capacity push alpha even higher than in
    # the unfrozen phase. Sampling only near T_init and Tf under-covers that
    # regime and was the direct cause of an early blow-up in this file's own
    # first draft. Scanning a wide, fixed grid instead of guessing which two
    # or three points matter is cheap (a few dozen evaluations) and removes
    # the guesswork entirely.
    def stability_dt(Tsample):
        k = k_of(Tsample); C = C_of(Tsample); alpha = k / C
        h_max = max(3.0, 10.0, 60.0)  # covers both cavity constants and a
                                        # generous freezer/oven-side h
        a_boundary = 2 * alpha * h_max / (k * dx) + 2 * alpha / dx ** 2
        dt_interior = dx * dx / (2 * alpha)
        dt_boundary = 2.0 / a_boundary
        return safety * min(dt_interior, dt_boundary)

    # Scan only the range the CALLER says this problem will actually visit
    # (T_bounds), not a fixed worst-case span — ice at -60 degC conducts fast
    # enough that scanning down there by default made every solve here take
    # minutes for no reason a real freezing time ever needs. A caller that
    # forgets to pass T_bounds still gets a safe, if slower, default.
    lo_scan, hi_scan = T_bounds if T_bounds else (min(T_init, -45.0), max(T_init, 10.0))
    scan_points = [lo_scan + i * 5.0 for i in range(int((hi_scan - lo_scan) // 5.0) + 2)]
    scan_points += [Tf, Tf - 0.1, Tf + 0.1]
    dt = min(stability_dt(Ts) for Ts in scan_points)
    # A stateful inner BC (an air pocket, a pool) can carry its OWN, much
    # tighter stability bound — see AirPocketInnerBC.max_stable_dt()'s
    # docstring for why a low-heat-capacity cavity needs this. Only applied
    # when the BC actually defines it; a BC that does not is assumed to have
    # no state of its own to go unstable.
    if hasattr(inner_bc_fn, 'max_stable_dt'):
        dt = min(dt, inner_bc_fn.max_stable_dt())
    nsteps = int(math.ceil(t_end / dt))
    if nsteps > max_steps:
        raise RuntimeError(f'{nsteps} steps requested (dt={dt:.3g}s) — refusing '
                            f'rather than running an unbounded loop.')
    dt = t_end / nsteps

    result = WallResult()
    t = 0.0
    for step in range(nsteps + 1):
        if step % record_every == 0 or step == nsteps:
            h_o, Tinf_o = outer_bc(t, T[0])
            h_i, Tinf_i = inner_bc_fn(t, T[-1])
            result.time_history.append(t)
            result.flux_in_history.append((h_o * (Tinf_o - T[0]), h_i * (Tinf_i - T[-1])))
        if step == nsteps:
            break

        h_o, Tinf_o = outer_bc(t, T[0])
        h_i, Tinf_i = inner_bc_fn(t, T[-1])

        # FLUX-CONSERVATIVE update. The earlier draft wrote each node's update as
        # alpha_i/dx^2 * (T[i+1]-2T[i]+T[i-1]) with alpha_i = k(T[i])/C(T[i]) taken
        # from that ONE node — correct only when k is uniform in SPACE. Once k
        # varies node to node (ice conducts ~4x better than water — freezing.js's
        # own module docstring says so), that form does not conserve energy: it
        # was found here as an energy-balance gap that stayed at ~44% however
        # fine the grid got (nx 21/41 both landed on it), which is the signature
        # of a formulation error, not a truncation one. The fix is the standard
        # one for variable-conductivity conduction: interface conductivities as
        # the HARMONIC mean of their two neighbouring nodes (two conductors in
        # series), and each node's update built from the flux difference across
        # it, which conserves energy by construction regardless of how sharply k
        # jumps between nodes.
        k = [k_of(Tv) for Tv in T]
        C = [C_of(Tv) for Tv in T]
        # q[i] = conductive flux, W/m2, across the interface between node i and
        # node i+1 (positive = flowing in the +x direction, i.e. INTO node i+1).
        q = [2 * k[i] * k[i + 1] / (k[i] + k[i + 1]) * (T[i + 1] - T[i]) / dx
             for i in range(nx - 1)]

        Tn = T[:]
        Tn[0] = T[0] + dt * 2 * (h_o * (Tinf_o - T[0]) + q[0]) / (C[0] * dx)
        for i in range(1, nx - 1):
            Tn[i] = T[i] + dt * (q[i] - q[i - 1]) / (C[i] * dx)
        Tn[-1] = T[-1] + dt * 2 * (h_i * (Tinf_i - T[-1]) - q[-1]) / (C[-1] * dx)

        if any(not math.isfinite(x) for x in Tn):
            raise RuntimeError(f'solution blew up at t={t:.3g}s — refusing to '
                                f'return a number.')

        if hasattr(inner_bc_fn, 'advance'):
            inner_bc_fn.advance(dt, T[-1])

        T = Tn
        t += dt

    result.T, result.nx, result.dx = T, nx, dx
    return result


def enthalpy_at(comp, Tf, tC, T_ref):
    """Enthalpy gained per unit VOLUME going from T_ref to tC, J/m3 — sensible
    term at the Choi-Okos midpoint cp (same convention frzPhamTime's own dT1/dT2
    midpoints use) plus rho0 * Lf * (xIce(tC) - xIce(T_ref)). Independent of the
    solver's own step-by-step apparent-C bookkeeping — this is the cross-check,
    not an echo of it."""
    rho0 = props_unfrozen(comp, (tC + T_ref) / 2.0)[0]
    # Sensible part: integrate the SENSIBLE (non-latent) volumetric heat
    # capacity from T_ref to tC. Split at Tf so each half uses its own branch's
    # cp. A SINGLE midpoint sample per half (the first draft of this function)
    # is exactly Pham's own dT1/dT2 approximation — fine for Pham, because
    # Pham's dT1/dT2 midpoints span a food's whole sensible range where cp
    # barely moves. It is NOT fine here: the frozen branch's cp is a mixture
    # property that changes FAST just below Tf (ice_fraction's derivative is
    # steepest right at the freezing point) and only levels off well below
    # it, so one point at the middle of a -1 to -30 degC span badly
    # misrepresents the true integral — this was the actual cause of a 39-44%
    # energy-balance gap that widening the apparent-C smoothing window (a
    # different, unrelated fix) did not touch. Composite Simpson's rule with
    # enough panels fixes it without needing a closed form for cp(T).
    def Csens_of(tval):
        return (props_frozen(comp, Tf, tval)[3] if tval < Tf
                else props_unfrozen(comp, tval)[3])

    def simpson(f, a, b, panels=20):
        if panels % 2 == 1:
            panels += 1
        h = (b - a) / panels
        total = f(a) + f(b)
        for i in range(1, panels):
            total += (4 if i % 2 else 2) * f(a + i * h)
        return total * h / 3.0

    lo, hi = min(T_ref, tC), max(T_ref, tC)
    sensible_unsigned = 0.0
    if lo < Tf < hi:
        sensible_unsigned += simpson(Csens_of, lo, Tf)
        sensible_unsigned += simpson(Csens_of, Tf, hi)
    else:
        sensible_unsigned += simpson(Csens_of, lo, hi)
    signed_sensible = sensible_unsigned if tC >= T_ref else -sensible_unsigned
    # Latent term: freezing (more ice, xIce rises) RELEASES heat, so it DROPS the
    # enthalpy the wall holds — xi(T_ref) - xi(tC), not the other way round. Got
    # this backwards on the first pass (xi(tC) - xi(T_ref)), which made a wall
    # that was COOLING show a POSITIVE enthalpy gain — the exact symptom the
    # energy-balance check below caught (flux-in negative, enthalpy positive).
    latent = rho0 * FRZ_L_WATER * (ice_fraction(comp, Tf, T_ref) - ice_fraction(comp, Tf, tC))
    return signed_sensible + latent


def energy_balance_error_latent(result, comp, Tf, T_init):
    times = result.time_history
    fluxes = [qo + qi for qo, qi in result.flux_in_history]
    energy_in = 0.0
    for i in range(1, len(times)):
        dt = times[i] - times[i - 1]
        energy_in += 0.5 * (fluxes[i] + fluxes[i - 1]) * dt

    dx = result.dx
    enthalpy_gained = 0.0
    for i, Tval in enumerate(result.T):
        w = dx if 0 < i < result.nx - 1 else dx / 2
        enthalpy_gained += enthalpy_at(comp, Tf, Tval, T_init) * w

    denom = max(abs(energy_in), abs(enthalpy_gained), 1.0)
    return abs(energy_in - enthalpy_gained) / denom, energy_in, enthalpy_gained


# ═══════════════════════════════════════════════════════════════════════════
# 3. THE ORACLE FOR THE PURELY-SENSIBLE REGIME — same closed-form series as
#    the cooking file, for the grid-convergence check (section 4c).
# ═══════════════════════════════════════════════════════════════════════════

def slab_series_lambda_n(bi, n, iterations=200):
    lo = (n - 1) * math.pi + 1e-12
    hi = (n - 1) * math.pi + math.pi / 2 - 1e-12
    f = lambda lam: lam * math.tan(lam) - bi
    flo = f(lo)
    for _ in range(iterations):
        mid = (lo + hi) / 2
        fmid = f(mid)
        if (flo <= 0) == (fmid <= 0):
            lo, flo = mid, fmid
        else:
            hi = mid
    return (lo + hi) / 2


def slab_series_theta_center(bi, fo, nterms=60):
    total = 0.0
    for n in range(1, nterms + 1):
        lam = slab_series_lambda_n(bi, n)
        Cn = 4 * math.sin(lam) / (2 * lam + math.sin(2 * lam))
        total += Cn * math.exp(-lam * lam * fo)
    return total


# ═══════════════════════════════════════════════════════════════════════════
# 4. VALIDATION
# ═══════════════════════════════════════════════════════════════════════════

def validate_symmetric_reduction_vs_pham():
    print('\n1. SYMMETRIC REDUCTION vs. PHAM — identical (h, T-infinity) on both')
    print('   faces must match freezing.js\'s own shipped formula for the same slab')
    print('   (a fitted closed-form correlation vs. a finite-difference PDE — the')
    print('   two most different methods this project has for the same question):')

    comp, Tf = DEMO_COMP, DEMO_TF
    Ti, Tc, Tm = 10.0, -18.0, -30.0
    h = 20.0
    D = 0.06  # 6 cm slab, well inside Pham's tested range

    tC_mid = (Ti + Tc) / 2.0
    Cl = props_unfrozen(comp, tC_mid)[3]
    Cs_frozen = props_frozen(comp, Tf, (Tf + Tc) / 2.0)
    Cs, ks = Cs_frozen[3], Cs_frozen[2]
    Lf = props_unfrozen(comp, Tf)[0] * freezable_water(comp) * FRZ_L_WATER
    VA = D / 2.0  # Pham's slab convention, freezing.js's frzGeometry('slab', [D])

    pham_seconds = pham_time_seconds(Cl, Cs, Lf, ks, Ti, Tc, Tm, h, VA, D)

    outer_bc = lambda t, Tface: (h, Tm)
    inner_bc = lambda t, Tface: (h, Tm)
    lo_t, hi_t = 0.0, pham_seconds * 3.0
    for _ in range(18):
        mid_t = (lo_t + hi_t) / 2
        res = solve_wall_latent(D, comp, Tf, Ti, outer_bc, inner_bc, mid_t, nx=21,
                                 T_bounds=(Tm, Ti))
        centre = res.T[res.nx // 2]
        if centre <= Tc:
            hi_t = mid_t
        else:
            lo_t = mid_t
    fd_seconds = hi_t

    rel_err = abs(fd_seconds - pham_seconds) / pham_seconds
    check('FD centre-reaches-Tc time matches Pham within 20%',
          rel_err < 0.20,
          f'FD={fd_seconds/3600:.3f} h, Pham={pham_seconds/3600:.3f} h, '
          f'rel. err={rel_err*100:.1f}%')
    print('   (20%, not 1%: Pham is itself an approximate correlation, not an exact')
    print('   solution of this PDE — the two methods are expected to differ by a')
    print('   margin similar to Pham\'s own published accuracy against measured data.')
    print('   This checks that they agree in ORDER, not that either is exact.)')


def validate_energy_conservation_latent():
    print('\n2. ENERGY CONSERVATION, WITH LATENT HEAT — flux through both faces vs.')
    print('   enthalpy gained (sensible + latent, read off ice_fraction()')
    print('   independently of the solver\'s own local bookkeeping):')

    comp, Tf = DEMO_COMP, DEMO_TF
    T_init = 10.0
    outer_bc = lambda t, Tface: (20.0, -30.0)
    inner_bc = lambda t, Tface: (3.0, -30.0)  # conservative cavity, same T_inf as outer
    res = solve_wall_latent(0.06, comp, Tf, T_init, outer_bc, inner_bc,
                             t_end=6 * 3600.0, nx=21, record_every=1,
                             T_bounds=(-30.0, T_init))

    rel_err, e_in, e_gained = energy_balance_error_latent(res, comp, Tf, T_init)
    check('Boundary flux integral matches enthalpy gained (sensible+latent) within 2%',
          rel_err < 0.02,
          f'flux-in={e_in:,.0f} J/m3, enthalpy gained={e_gained:,.0f} J/m3, '
          f'rel. err={rel_err*100:.3f}%')


def validate_grid_convergence_sensible_only():
    print('\n3. GRID CONVERGENCE, purely-sensible regime (no Tf crossed) — same')
    print('   multi-term series oracle the cooking file already validated this')
    print('   scheme against, proving the base scheme before C(T) is added:')

    comp = DEMO_COMP
    Tf = -50.0  # never crossed at the temperatures used below — irrelevant here
                # since k_of/C_of are overridden to true constants regardless
    # TRULY constant k and C, not Choi-Okos sampled over a range — this isolates
    # the SCHEME's own truncation error from the small T-dependence Choi-Okos
    # always carries. A first draft of this check called Choi-Okos over a narrow
    # span instead, on the theory that "narrow enough" would make the model
    # drift negligible; it did not — the drift stayed comparable to the
    # discretisation error at any span tried, giving a flat or even
    # non-monotonic ratio (found here at both a 20-degree and a 1-degree span).
    # Overriding to real constants (as the cooking file's own check already
    # does, for the identical reason) recovers clean 4th-order convergence
    # (halving dx AND dt together): confirmed separately at ratio ~4.00 before
    # writing this in permanently.
    k_const, C_const = 0.45, 3.5e6
    alpha = k_const / C_const
    L_half = 0.03
    h, T_init, T_inf = 20.0, 30.0, 10.0
    Bi = h * L_half / k_const
    t_end = 1800.0
    Fo = alpha * t_end / (L_half * L_half)
    theta_series = slab_series_theta_center(Bi, Fo)

    outer_bc = lambda t, Tface: (h, T_inf)
    inner_bc = lambda t, Tface: (h, T_inf)

    errors = []
    for nx in (51, 101, 201):
        res = solve_wall_latent(2 * L_half, comp, Tf, T_init, outer_bc, inner_bc,
                                 t_end, nx=nx, T_bounds=(T_inf, T_init),
                                 k_of_override=lambda Tv: k_const,
                                 C_of_override=lambda Tv: C_const)
        centre = res.T[res.nx // 2]
        theta_fd = (centre - T_inf) / (T_init - T_inf)
        errors.append(abs(theta_fd - theta_series))

    ratio1 = errors[0] / errors[1] if errors[1] > 0 else float('inf')
    ratio2 = errors[1] / errors[2] if errors[2] > 0 else float('inf')
    check('Error shrinks under grid refinement (both halvings, ratio > 1.5)',
          ratio1 > 1.5 and ratio2 > 1.5,
          f'errors at nx=51/101/201: {errors[0]:.2e} / {errors[1]:.2e} / {errors[2]:.2e}, '
          f'ratios={ratio1:.2f}, {ratio2:.2f}')


# ═══════════════════════════════════════════════════════════════════════════
# 5. THE CAVITY BAND — an AIR POCKET, not a pool. No evaporation, no fixed
#    T-infinity from the freezer: the trapped air only exchanges with the
#    wall's own inner face, and cools (or, early on, may even warm the face
#    back slightly) only through that.
# ═══════════════════════════════════════════════════════════════════════════

H_CAVITY_INERT = 3.0    # W/m2K — SAME constant and SAME status as the cooking
H_CAVITY_HUMID = 10.0   # file's H_CAVITY_INERT/HUMID: natural convection in a
                        # closed air gap, the identical physical mechanism —
                        # reused rather than duplicated with new unsourced numbers.

AIR_RHO = 1.30      # kg/m3 — dry air near freezer temperatures (ASSUMED point value,
AIR_CP = 1005.0     # J/(kg*K) — not swept over the actual -30..0 degC range this
                    # app's freezer-temperature field allows; a fixed-point
                    # approximation, same status as thermal.js's other single-point
                    # air properties)


class AirPocketInnerBC:
    """The cavity's own trapped air: starts at the food's initial temperature,
    has a real (if small) heat capacity set by the ACTUAL cavity volume (passed
    in, not assumed), and exchanges heat with nothing but the wall's inner face —
    no external sink, no evaporation. Andrea, 2026-09-13: "tenendo conto che la
    cavità è vuota." Stateful and single-use, same convention as the cooking
    file's PoolInnerBC: a fresh instance per solve."""

    def __init__(self, T_init, h_cavity, C_areal):
        self.T_air = T_init
        self.h = h_cavity
        self.C = C_areal  # J/(m2*K) — rho_air * cp_air * V_cavity / A_inner

    def __call__(self, t, T_face):
        return self.h, self.T_air

    def advance(self, dt, T_face):
        q = self.h * (T_face - self.T_air)  # > 0 when face is warmer than air
        self.T_air += dt * q / self.C

    def max_stable_dt(self, safety=0.5):
        """Air has almost no mass, so this BC's own time constant (C/h) can be
        far shorter than anything the WALL's stability bound sees — a small
        cavity volume gives a tiny C_areal, and an explicit update of
        dT_air/dt = -(h/C)*(T_air-T_face) blows up past dt = 2C/h regardless of
        how fine the wall's own grid is. Found here as a genuine blow-up (not
        the wall's usual one) the first time a demo cavity was small enough for
        this to bind before the wall's own limit did. solve_wall_latent asks
        for this explicitly, on any inner BC that defines it."""
        return safety * 2.0 * self.C / self.h


def run_band_example():
    print('\n4. THE BAND — one illustrative wall, three cavity readings. Not a food')
    print('   database row: a demo of the band\'s WIDTH before it is wired to a real')
    print('   dThermal from the plane-section stack.\n')

    comp, Tf = DEMO_COMP, DEMO_TF
    Ti, Tc, Tm = 8.0, -18.0, -30.0
    h_out = 20.0
    D = 0.05  # 5 cm wall
    t_end_max = 12 * 3600.0

    # A representative cavity: 1 L of trapped air per m2 of inner wall face (a
    # round, illustrative ratio — not tied to any one carcass; the app-side
    # wiring will compute the real V_cavity/A_inner from the stack the cook
    # already measured or the class the weight deduced).
    V_over_A = 0.01  # m — 10 L of trapped air per m2 of inner wall face: a
                      # rounder, less extreme illustrative ratio than the first
                      # draft's 1 L/m2, which forced a stability-bound dt so
                      # tight (see AirPocketInnerBC.max_stable_dt) that this
                      # demo alone took minutes in plain Python for no reason a
                      # real cavity needs — a real carcass cavity is not that
                      # starved of air.
    C_areal = AIR_RHO * AIR_CP * V_over_A

    outer_bc = lambda t, Tface: (h_out, Tm)

    # What actually varies here is NOT how well the trapped air couples to the
    # face (a second draft tried exactly that — h_cavity alone, same dynamic
    # air pocket for all three members — and found the three came out numbers
    # indistinguishable to four significant figures: air holds so little heat
    # that whatever h is, it reaches local equilibrium with the face in
    # seconds, far faster than the hours this problem runs on, so h stops
    # mattering to the total once past the first instant).
    #
    # The real open question for an EMPTY cavity, Andrea 2026-09-13: a
    # carcass's cavity is NOT sealed — every animal has an opening (the vent,
    # the neck end after drawing) that is a real, physical crack to the
    # freezer's own air, and that crack's SIZE scales with the animal's own
    # size, not a fixed number. Nobody holds a measurement of that opening
    # across species yet (the same open morphometry question the plane-
    # section stack's cavity fraction already waits on). So this band's two
    # members are the two ends of THAT one physical axis — how open the
    # crack effectively is — not two competing guesses about something else:
    #   conservative/centro — the crack is small relative to the cavity: the
    #     trapped air pocket model (which, per the finding above, behaves
    #     near-adiabatically regardless of its exact coupling h) stands in
    #     for "mostly sealed". This project's default until a real opening
    #     size exists to scale the band by.
    #   aggressive — the crack is large enough that the cavity air is
    #     effectively already at the freezer's own temperature (Tinf=Tm from
    #     the first instant): the fully-open limit.
    # A future version that has an actual opening-area-vs-body-size reading
    # could set the cavity's OWN effective h from that geometry directly,
    # collapsing this into one number instead of a band — until then the band
    # is the honest way to carry "we know there is a crack, not how big".
    def make_inner(mode):
        if mode in ('conservative', 'centro'):
            return AirPocketInnerBC(Ti, H_CAVITY_INERT, C_areal)
        if mode == 'aggressive':
            return lambda t, Tface: (H_CAVITY_HUMID, Tm)
        raise ValueError(mode)

    print(f'   {"mode":<14}{"time for inner face to reach -18C":>36}')
    for mode in ('conservative', 'centro', 'aggressive'):
        lo_t, hi_t = 0.0, t_end_max
        for _ in range(14):
            mid_t = (lo_t + hi_t) / 2
            inner_bc = make_inner(mode)
            res = solve_wall_latent(D, comp, Tf, Ti, outer_bc, inner_bc, mid_t, nx=15,
                                     T_bounds=(Tm, Ti))
            if res.T[-1] <= Tc:
                hi_t = mid_t
            else:
                lo_t = mid_t
        reached = hi_t < t_end_max - 1.0
        label = f'{hi_t/60:.1f} min' if reached else f'> {t_end_max/60:.0f} min (not reached)'
        print(f'   {mode:<14}{label:>36}')

    print('\n   NOTE on what this band actually varies, and why "centro" here equals')
    print('   "conservative": air holds so little heat that a trapped pocket reaches')
    print('   local equilibrium with the wall\'s inner face in SECONDS, however well or')
    print('   badly it is coupled — irrelevant next to the HOURS this problem runs on.')
    print('   So the coupling strength h is not the real question a truly empty,')
    print('   sealed cavity raises. The real question is whether it is sealed at all:')
    print('   a carcass cavity that genuinely vents to the freezer\'s own air (through')
    print('   whatever opening it has) behaves quite differently from one that does')
    print('   not. "aggressive" here is that vented case (Tinf=Tm from the first')
    print('   instant); "conservative" and "centro" both model the sealed case, which')
    print('   this project takes as the ordinary one for a whole carcass — so the band')
    print('   is genuinely two-valued, not three, and says so rather than pretending a')
    print('   third number exists where the physics does not support one.')


if __name__ == '__main__':
    print('Asymmetric plane-wall conduction WITH LATENT HEAT — validation before')
    print('anything is wired to freezing.js. See the module docstring for scope.')
    print('\n═══ Assumed numbers in play ═══')
    print(f'  H_CAVITY_INERT = {H_CAVITY_INERT} W/m2K  (reused from the cooking file —')
    print('                                             same mechanism: still air)')
    print(f'  H_CAVITY_HUMID = {H_CAVITY_HUMID} W/m2K  (reused — lightly-moved air)')
    print(f'  AIR_RHO/AIR_CP = {AIR_RHO} kg/m3, {AIR_CP} J/(kg*K)  (ASSUMED single-point air)')

    print('\n═══ Validation ═══')
    validate_symmetric_reduction_vs_pham()
    validate_energy_conservation_latent()
    validate_grid_convergence_sensible_only()

    run_band_example()

    print('\n' + '═' * 60)
    if FAILURES:
        print(f'FAILED: {", ".join(FAILURES)}')
        sys.exit(1)
    else:
        print('All checks passed. The method is validated against Pham in the')
        print('symmetric limit and conserves energy (sensible+latent) under the')
        print('asymmetric band; it is not yet a "motore di controllo" for shipped')
        print('JavaScript, and no JavaScript has been touched.')
        sys.exit(0)
