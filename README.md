# Plane sections and cavities: a new method for an old kitchen problem

<blockquote style="font-family: Georgia, 'Times New Roman', serif; font-size: 0.95em; font-style: italic; border-left: 3px solid #999; padding-left: 1em; margin-left: 0; color: #444;">
«Philosophy is written in this grand book, the universe, which stands
continually open to our gaze. But the book cannot be understood unless one
first learns to comprehend the language and read the letters in which it is
composed. It is written in the language of mathematics, and its characters
are triangles, circles, and other geometric figures, without which it is
humanly impossible to understand a single word of it; without these, one
wanders about in a dark labyrinth.»
<br>— Galileo Galilei, <em>The Assayer</em> (1623), translated by Stillman Drake.
</blockquote>

**Why this idea is public.** Andrea Petrucci, the method's originator:

> "We publish it: we think we are the first to do so, since we have found no
> scientific literature on it, and knowledge sets us free and is a common
> good."

*Sapere Aude* — dare to know, dare to use your own reason. **No one should be
able to profit privately from an idea that should stay public**: science
exists to be shared, not locked in a patent.

**Every figure in this document is checkable by anyone** — no code is shown
here, only the mathematics and the results, reproducible with the scripts
named at the end of each section.

---

## 1 · The problem, in one sentence

To estimate how long a piece of food takes to cook or freeze, the software
first has to know **what shape it has**. For an irregular piece — a whole
chicken, a beef half, a whole fish — the only practical choice so far was to
close it inside an **infinite slab** (safe, but a considerable overestimate)
or an imaginary **box**. On an 8 cm cylinder, the box around it holds **27.3%
more** volume than the true one — not a measurement error, plain arithmetic:
a box around any round body over-states it by exactly 4/π − 1, and around a
sphere by **91%**.

## 2 · The method: plane sections, not a 3D model

Instead of one shape chosen in advance, the piece is imagined sliced into
thin, regular sections — like a sliced salami — each measured (or, where it
cannot be, cautiously inferred from weight). The volume follows from summing
them by **Cavalieri's principle**: the same integral used for a cone or a
sphere, applied slice by slice to the real shape instead of an assumed one.

On the same 603.19 cm³ cylinder (true value, exact formula):

| Method | Computed volume | Error from the true value |
|---|---|---|
| Box | 768.00 cm³ | **+27.3%** |
| Plane sections | 602.80 cm³ | **−0.06%** |

Over **450 times more precise**, on the identical reading. Confirmed too on
an oblate ellipsoid (1340.07 vs. 1340.41 cm³, −0.03%) and a hollow tube (wall
thickness 2.99949 vs. 3.00000 cm, −0.02%).

*Reproducible: `node tests/plane_sections_report.js`, in the application's repository (private today — not yet extracted standalone into this public repository).*

## 3 · What changes for a real FREEZING time

Comparing the **infinite slab** — the geometry the scientific community has
used and approved for decades for an irregular piece, via Pham's (1986)
method — against plane sections, on the identical readings:

| Piece | Slab (standard method) | Plane sections (our method) | Difference |
|---|---|---|---|
| 8 cm cylinder | 7.09 h | 3.02 h | **−57%** |
| Tied roast | 11.86 h | 4.76 h | **−60%** |
| Beef half | 20.57 h | 8.45 h | **−59%** |
| Whole fish | 10.59 h | 3.83 h | **−64%** |

The estimated time **drops by 57% to 64%** against the standard method.
Important: this shows the **geometry** is more faithful — it does not, by
itself, prove the **time** is too. See §7.

## 4 · The new problem: the HOLLOW carcass

A whole chicken, turkey or rabbit is not solid: it has an internal cavity.
Pham's equations were fitted on **solid** bodies — using them on a hollow one
would mean applying a formula outside its own domain, and this app refuses to
do that silently. A different method was needed for two distinct problems —
**cooking** (the oven heats) and **freezing** (the freezer cools) — both
solved with the same principle.

### The shared principle: ASYMMETRIC heating/cooling

**From the measurements to the temperatures**, the path is: measure (or
deduce from weight) the carcass's sections → derive the meat thickness
between the outer surface and the cavity → solve heat conduction through that
thickness, **with two different conditions on the two faces of the same
wall** — the oven or the freezer on one side, the internal cavity on the
other.

Classical models (Cleland, Pham) assume the **same** temperature and the
same heat-transfer coefficient over the whole surface of the piece — a
convenient assumption that works for a solid piece, but which a cavity breaks
on purpose. Our engine solves a different physical problem: conduction with
**asymmetric** boundary conditions, verified with a numerical method (finite
differences) independent of the classical formula.

- **In cooking**: the cavity heats up and — where the carcass allows it — a
  pool of water forms inside it, evaporating as it warms (never fixed 100°C
  steam from the first instant: a real energy balance, corrected twice live
  by Andrea before reaching this version).
- **In freezing**: the cavity is **empty** (air, not a pool) — no
  evaporation. The real uncertainty is not how well the air couples to the
  wall (irrelevant: air holds so little mass it equilibrates in seconds,
  against the hours the problem runs on), but whether the natural opening
  every animal has (which scales with its own size) lets the freezer's own
  air in or not.

## 5 · The new freezing engine — validated three independent ways

| Check | What it verifies | Result |
|---|---|---|
| **Against Pham** (symmetric case) | With identical conditions on both faces, the new method must return the same time as the standard method already shipped | **9.5% difference** (threshold 20%: these are two radically different methods — a closed-form correlation against a numerical simulation) |
| **Energy conservation** | Heat entering and leaving, integrated over time, must match the energy actually gained or lost (sensible heat + the latent heat of the ice forming) | **1.49% difference** (threshold 2%) |
| **Grid convergence** | Refining the computational grid should shrink the error predictably | **Exactly 4x smaller** at each refinement — the theoretically perfect rate |

*Reproducible: `python asymmetric_wall_freezing_conduction.py`, included in this same repository (standalone script, standard library only).*

## 6 · Two examples, identical measurements — one solid carcass, one hollow

A **6 cm** wall in both cases, same starting temperature (4°C), same freezer
(−30°C), same target (−18°C):

| Carcass | Time to reach −18°C |
|---|---|
| **Solid** (Pham, standard method) | **6.15 h** |
| **Hollow, sealed opening** (our prudent case) | **12.61 h** |
| **Hollow, with a real crack to the freezer** | **7.17 h** |

**Why the hollow one takes longer, at the same wall thickness.** In a solid
piece, cold arrives **from both sides** toward the centre, effectively
halving the path. In a hollow piece with a sealed cavity, cold arrives **from
one side only**: the inner face has to wait for the whole thickness to cool,
with no help from the other side. If the carcass has a real crack that lets
cold air in, the time moves closer to that of a solid piece. **This is
today's single most important open question**: how large is that crack, in a
real animal?

## 7 · The honesty that matters most

**Knowing the geometry is more precise does not mean knowing the time is
too.** The formulas that turn a shape into a time were fitted on idealised
solids — giving them a truer shape does not automatically guarantee a truer
time. This app has already been bitten by exactly this: **in August, a more
realistic pair of shape factors made freezing times up to 31% too short**
before the mistake was found — that time, in the dangerous direction.

So: the geometry is **measurably** better (§2), checkable with a calculator.
Whether the time is too, **we do not know yet** — that needs a probe in a
real roast or a real frozen carcass, and we have not done that verification.

## 8 · The sources consulted — including the discarded ones

To build the cavity model we wrote, on **11 September 2026**, to the
Bertoliana Civic Library in Vicenza, Italy, asking for **13 sources** —
articles and textbooks. Outcome, reported honestly:

- **Found and freely consultable, in person**: Singh & Heldman, *Introduction
  to Food Engineering* (La Vigna Library, Vicenza).
- **Found but only for a fee** (inter-library loan, €10): Datta, *Biological
  and Bioenvironmental Heat and Mass Transfer* — **discarded**, on a stated
  principle: free sources only.
- **Unavailable on principle**: the complete ASHRAE Handbook – Refrigeration,
  purchasable only directly from the publisher — **discarded**.
- **Unavailable through the circuit**: Polley, Snyder & Kotnour (1980),
  consultable only at the International Institute of Refrigeration's own
  library — **discarded**.
- **No reply from the library** on 7 further items (Cleland 1987, Cleland &
  Earle from the 1980s, Hossain et al. 1992, Siripon/Tansakul/Mittal 2007,
  Incropera & DeWitt, Sun, Rahman).

**In parallel**, three direct requests to the authors themselves:
- **G. S. Mittal** (30 July 2026) and **Ampawan Tansakul** (28 August 2026),
  co-authors of the only published precedent that treats a whole chicken as a
  stack of cross-sections (Siripon, Tansakul & Mittal, 2007) — **no reply
  from either**.
- **Ors Petnehazy** (11 September 2026), author of the only public CT
  tomography of a turkey's cavity — email sent, delivery unconfirmed. His
  **2022 dataset (CC BY, free) has already been downloaded and used**, no
  permission required: it is the only real source of animal morphometry
  actually in use in this work today.

**Honest total: 16 sources sought for this specific piece of work (13
through the library + 3 direct to authors), 2 genuinely available for free, 1
public dataset in use, everything else absent or discarded.** A real
carcass's cavity remains, for now, a declared assumption — not a
measurement.

## 9 · Control engines — what they are, how many we have

**A control engine does not repeat a calculation: it checks it by a
different route.** A second calculation by the same method is not a check,
it is an echo. This document's own examples (§5) are the textbook case: the
new freezing engine was checked against (a) a completely different formula
(Pham), (b) a conservation law that knows nothing about our model (energy),
(c) the known behaviour of a numerical scheme under grid refinement.

**How many we have**: the whole app is covered by **82 automated test
suites, 20,150 assertions**, re-run on every change and mechanically checked
against the figures the project's own documents claim (a check on the
checks). Among these, a core of over **26 genuine control engines** — each
catalogued with its own history, often born from a real bug found and fixed
— covers geometry, thermal transfer, freezing, lethality/safety, and the
app's own claims about itself.

*Reproducible: `node tests/run_all.js` from the application repository's root (private today).*

## 10 · Why we are publishing now

**As far as we know, no public source applies the repeated-plane-sections
method — with a deduced or measured cavity, for both heating and cooling —
to a cooking or freezing time in a kitchen.** We are not claiming credit for
an old idea (Cavalieri, 17th century): we are asking that it be on the
record, with a firm date, that we were first to write it down and publish it
applied to this specific problem. No patent: knowledge is a common good, and
no one should be able to make a private profit from it.

## How to cite this work

KitchenMetrics is **free**, with no advertising: there is no commercial
motive behind this request. If you use the repeated-plane-sections-with-
cavity method (deduced or measured) — or this text — in your own work, a
project, a paper, another app, crediting **Andrea Petrucci** as the
originator would simply be **a matter of intellectual honesty**: not an
obligation we can enforce on a mathematical idea (copyright protects
expression, not method — see `LICENSE-SCOPE.md`), just what one would expect
between people acting in good faith.

> Petrucci, A. (2026). *Repeated plane sections with a deduced or measured
> cavity, applied to cooking and freezing time.* KitchenMetrics.
> https://github.com/apetrucci90-rgb/KitchenMetrics

---

*Published 13 September 2026. Every figure is reproducible from the scripts
named in this public repository. Italian version:
`POST-GITHUB-SEZIONI-PIANE-IT.md`.*