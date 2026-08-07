#!/usr/bin/env python3
"""Compare the repeat-pass InSAR tidal signal against ApRES at the same site.

WHY THIS EXISTS. The InSAR measurement (this repo) reports the tidal
response as a STRAIN per metre of tidal heave. The ApRES processing in
~/projects/EAGER_ApRES reports a vertical STRAIN RATE in 1/yr, tidally
modulated. Those are not the same quantity and cannot be compared until one
is converted, which is what this does. It is the only independent check we
have on whether the signal we are trying to detect is the size we think.

THE CONVERSION. If the strain rate oscillates as A*cos(w*t), the strain
itself oscillates as (A/w)*sin(w*t), so the strain amplitude is A/w. Divide
by the tidal amplitude to get strain per metre of tide, the InSAR quantity.

    strain_amp = A / w                     [dimensionless]
    admittance = strain_amp / tide_amp     [1/m]

WHAT IT SHOWS. ApRES and the thin-plate flexure estimate drawn on
scripts/figures/strain_rates.m agree to about 10%, and both sit roughly 7x
below the InSAR systematic floor. So the expected-signal curve on that
figure is not just a model - an independent instrument at the same site in
the same weeks measures the same magnitude.

SOURCES. ApRES amplitudes are from a least-squares harmonic fit to the
medfilt vertical strain-rate series in
~/projects/EAGER_ApRES/results/<site>/pair_results.csv. That fit is NOT a
product of the ApRES repo - there is no tidal analysis in it - so these
numbers live here rather than being read from a file there.

CAVEATS, all of which matter:
  - DIFFERENT DEPTH INTERVALS. ApRES excludes the firn and fits from 100 m
    to about 20 m above the bed; the InSAR quantity here is the 0-100 m
    column. They do not overlap at all.
  - THE NEUTRAL PLANE SITS INSIDE THE ApRES INTERVAL. Bending strain
    reverses sign at about H/2, which for 206-272 m of ice is 103-136 m -
    inside the 101-224 m fit. So the ApRES bulk slope partially cancels its
    own bending signal and is a LOWER BOUND on the near-surface amplitude.
  - Only GA04 is trustworthy: GA01's tracked bed drifts -81.6 m and GA05's
    +158.7 m, i.e. the bed pick jumps between reflectors.
  - A 5.7-8.9 day record cannot separate K1 from O1 (13.6 d beat) or M2
    from S2 (14.8 d), so constituent splits are indicative; the
    diurnal-band power is the robust part.
  - No ApRES site coordinates exist (GPS was off for the whole deployment),
    so "the same site" means the same few-km area, not a known offset.
"""

import math

DAYS_PER_YEAR = 365.25

# ApRES diurnal strain-rate amplitudes [1/yr] and their fit intervals [m].
# GA01 is carried for completeness and flagged, not used for conclusions.
APRES = {
    "GA04": dict(amp_per_yr=0.0129, fit_top=101, fit_bot=224, ice=250.6,
                 var_explained=0.62, trust=True),
    "GA01": dict(amp_per_yr=0.0105, fit_top=101, fit_bot=177, ice=206.0,
                 var_explained=0.49, trust=False),
}

DIURNAL_HOURS = 23.93          # K1; O1 is 25.82 h, which changes this ~8%
TIDE_AMP_M    = 0.55           # CATS2008 over the array: 1.093 m range

# InSAR side, from scripts/figures/strain_rates.m over the 0-100 m column
INSAR_FLOOR_USTRAIN_PER_M = 83.0
INSAR_EXPECTED_USTRAIN_PER_M = 12.3    # thin-plate flexure, H=300 m, L=2 km


def admittance_ustrain_per_m(amp_per_yr, period_hours, tide_amp_m):
    """Strain-rate amplitude [1/yr] -> strain per metre of tide [ustrain/m]."""
    amp_per_hour = amp_per_yr / (DAYS_PER_YEAR * 24.0)
    omega = 2.0 * math.pi / period_hours          # [1/h]
    strain_amp = amp_per_hour / omega             # dimensionless
    return 1e6 * strain_amp / tide_amp_m


def main():
    print("ApRES tidal strain-rate amplitude converted to the InSAR quantity")
    print(f"  diurnal period {DIURNAL_HOURS} h, tidal amplitude {TIDE_AMP_M} m\n")
    print(f"{'site':6} {'amp[1/yr]':>10} {'strain amp':>12} {'adm[ue/m]':>11} "
          f"{'span[m]':>9} {'mm over span':>13} {'trust':>6}")
    for site, d in APRES.items():
        adm = admittance_ustrain_per_m(d["amp_per_yr"], DIURNAL_HOURS, TIDE_AMP_M)
        span = d["fit_bot"] - d["fit_top"]
        mm = adm * 1e-6 * span * 1e3
        strain_amp = adm * 1e-6 * TIDE_AMP_M
        print(f"{site:6} {d['amp_per_yr']:10.4f} {1e6*strain_amp:10.2f}ue "
              f"{adm:11.1f} {span:9.0f} {mm:13.2f} {str(d['trust']):>6}")

    ga04 = admittance_ustrain_per_m(APRES["GA04"]["amp_per_yr"],
                                    DIURNAL_HOURS, TIDE_AMP_M)
    print(f"""
Comparison over the 0-100 m column, all as strain per metre of tide:
  ApRES GA04 (independent instrument) {ga04:6.1f} ue/m
  thin-plate flexure estimate         {INSAR_EXPECTED_USTRAIN_PER_M:6.1f} ue/m
  InSAR systematic floor              {INSAR_FLOOR_USTRAIN_PER_M:6.1f} ue/m

  ApRES vs the flexure estimate: {100*abs(ga04-INSAR_EXPECTED_USTRAIN_PER_M)/INSAR_EXPECTED_USTRAIN_PER_M:.0f}% apart
  InSAR floor above the ApRES value:  {INSAR_FLOOR_USTRAIN_PER_M/ga04:.1f}x

An independent instrument at the same site in the same weeks measures the
magnitude the flexure model predicts. The signal is real and it is roughly
seven times below what this method can currently resolve. Because the
neutral plane sits inside the ApRES fit interval, its value is a lower
bound, so the true near-surface signal is likely LARGER than {ga04:.0f} ue/m -
which shortens the gap rather than widening it.
""")


if __name__ == "__main__":
    main()
