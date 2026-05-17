"""Apply pilot PASS / MARGINAL / FAIL rules from plan §2.

Inputs: per-system SystemStats keyed by `bindingsite_wt` and `bindingsite_pack`.
Output: a structured decision dict ready to dump to pilot_decision.json.

Decision rule:
  PASS     if R ≥ R_PASS AND tau_KM(WT) ≥ TAU_WT_PASS_NS AND f_ee(WT) ≤ F_EE_PASS_MAX
  FAIL     if R < R_MARGINAL OR f_ee(WT) > F_EE_MARGINAL_MAX OR tau_KM(WT) < TAU_WT_MIN_NS
  MARGINAL otherwise
"""
from __future__ import annotations
from dataclasses import asdict, dataclass

from .stats import SystemStats
from ..config import (
    F_EE_MARGINAL_MAX, F_EE_PASS_MAX, PACK_SYSTEM, R_MARGINAL, R_PASS,
    TAU_WT_MIN_NS, TAU_WT_PASS_NS, WT_SYSTEM,
)


@dataclass
class PilotDecision:
    call: str                # "PASS" | "MARGINAL" | "FAIL"
    R: float                 # tau_KM(WT) / tau_KM(pack); inf if pack-tau is 0
    tau_KM_wt: float
    tau_KM_pack: float
    f_ee_wt: float
    reasons: list[str]
    next_step: str

    def to_dict(self) -> dict:
        return asdict(self)


def _ratio(num: float, den: float) -> float:
    if den <= 0:
        return float("inf") if num > 0 else 0.0
    return num / den


def evaluate(stats_by_system: dict[str, SystemStats]) -> PilotDecision:
    if WT_SYSTEM not in stats_by_system or PACK_SYSTEM not in stats_by_system:
        raise KeyError(f"need both {WT_SYSTEM!r} and {PACK_SYSTEM!r} in stats_by_system")
    wt = stats_by_system[WT_SYSTEM]
    pk = stats_by_system[PACK_SYSTEM]
    R = _ratio(wt.tau_KM, pk.tau_KM)

    fail_reasons = []
    if R < R_MARGINAL:
        fail_reasons.append(f"R={R:.2f} < {R_MARGINAL} (no discrimination)")
    if wt.f_ee > F_EE_MARGINAL_MAX:
        fail_reasons.append(f"f_ee(WT)={wt.f_ee:.2f} > {F_EE_MARGINAL_MAX} (WT unstable in equilibration)")
    if wt.tau_KM < TAU_WT_MIN_NS:
        fail_reasons.append(f"tau_KM(WT)={wt.tau_KM:.2f} ns < {TAU_WT_MIN_NS} ns (WT does not bind)")

    if fail_reasons:
        return PilotDecision(
            call="FAIL", R=R, tau_KM_wt=wt.tau_KM, tau_KM_pack=pk.tau_KM,
            f_ee_wt=wt.f_ee, reasons=fail_reasons,
            next_step="Pivot to unbinding metaD or steered-MD with Jarzynski.",
        )

    pass_reasons = []
    if R >= R_PASS:
        pass_reasons.append(f"R={R:.2f} ≥ {R_PASS}")
    if wt.tau_KM >= TAU_WT_PASS_NS:
        pass_reasons.append(f"tau_KM(WT)={wt.tau_KM:.2f} ns ≥ {TAU_WT_PASS_NS} ns")
    if wt.f_ee <= F_EE_PASS_MAX:
        pass_reasons.append(f"f_ee(WT)={wt.f_ee:.2f} ≤ {F_EE_PASS_MAX}")

    if (R >= R_PASS and wt.tau_KM >= TAU_WT_PASS_NS and wt.f_ee <= F_EE_PASS_MAX):
        return PilotDecision(
            call="PASS", R=R, tau_KM_wt=wt.tau_KM, tau_KM_pack=pk.tau_KM,
            f_ee_wt=wt.f_ee, reasons=pass_reasons,
            next_step="Lock the protocol. Write the full 13-system calibration plan.",
        )

    marginal_reasons = []
    if R < R_PASS:
        marginal_reasons.append(f"R={R:.2f} in [{R_MARGINAL}, {R_PASS})")
    if wt.f_ee > F_EE_PASS_MAX:
        marginal_reasons.append(f"f_ee(WT)={wt.f_ee:.2f} in ({F_EE_PASS_MAX}, {F_EE_MARGINAL_MAX}]")
    if wt.tau_KM < TAU_WT_PASS_NS:
        marginal_reasons.append(f"tau_KM(WT)={wt.tau_KM:.2f} ns in [{TAU_WT_MIN_NS}, {TAU_WT_PASS_NS}) ns")

    return PilotDecision(
        call="MARGINAL", R=R, tau_KM_wt=wt.tau_KM, tau_KM_pack=pk.tau_KM,
        f_ee_wt=wt.f_ee, reasons=marginal_reasons,
        next_step=("Diagnose: re-run WT at ±25% force, extend t_max to 200 ns, "
                   "inspect 3 trajectories visually. Re-evaluate same criterion."),
    )
