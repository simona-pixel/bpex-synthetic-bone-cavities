# Distinct numerical tolerances

Inspection of the preserved source code during revision established two different historical numerical rules. The vertex-sampled guard uses unsigned distance greater than or equal to the requested clearance minus 0.0001 mm. The subsequently added continuous auditor uses continuous distance plus 0.000000001 mm greater than or equal to the requested clearance. Both core implementations are preserved.

The initial new wrapper's auxiliary fixed_0150_pass field used the sampled tolerance. That metadata field was corrected to the continuous auditor rule; the generator, distance values and final_criteria_pass were not affected. Aggregation recomputes the fixed-threshold decision from the stored continuous distance and does not depend on that auxiliary field. This correction does not replace any executed configuration or discard an adverse result.

Mathematical descriptions in this revision must state the two tolerances separately. Neither tolerance is a certified error bound of floating-point geometry.
