"""Compare retained F2 formal jets with two documented metadata corrections.

Only the old seven-coordinate gauge label and the copied GF(31) sentence are
normalized. Field values, prescribed coordinates, coefficients, residuals
and digests remain part of the exact comparison. No record is overwritten.
"""

REGULAR_GAUGE_VARIABLES = (
    "P_-5_d0",
    *tuple(f"P_7_d{degree}" for degree in range(2, 8)),
)


def boundary_text(prime: int) -> str:
    return (
        f"a finite formal jet over GF({prime}) does not specialize safely at "
        "lambda=1 and is not a finite-field point or characteristic-zero lift"
    )


def comparison_record(payload: dict) -> dict:
    """Return the mathematical record, retaining every arithmetic field."""
    record = {key: value for key, value in payload.items() if key != "software"}
    if record.get("claim_boundary") == boundary_text(31):
        record["claim_boundary"] = boundary_text(record["field"]["prime"])
    gauge = record.get("higher_order_gauge", {})
    if (
        gauge.get("name") == "seven-pole-coordinate-zero"
        and gauge.get("variables_prescribed_zero_from_order_two")
        == list(REGULAR_GAUGE_VARIABLES)
    ):
        record["higher_order_gauge"] = {**gauge, "name": "selected-coordinate-zero"}
    return record
