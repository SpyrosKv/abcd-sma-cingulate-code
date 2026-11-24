# Reporting conventions

## Primary and sensitivity networks

The primary analysis uses pooled residualisation and separate between-person networks at each visit. The wave-specific residualisation analysis is a sensitivity analysis using the corresponding complete-case participants. Its models do not fit a family random intercept with one participant per family at a visit. A retained family component in the pooled models effectively describes stable participant variation across repeated assessments.

Primary full-wave networks and paired comparison networks use different participant sets. An edge weight in a paired comparison must not be described as the corresponding full-wave network estimate.

## Permutations and corrections

The six revised primary paired comparisons use 10,000 permutations. Their raw Monte Carlo p-values use (exceedances + 1) / (permutations + 1), and BH/Holm corrections apply within each family of 153 edges. These p-values are estimates with finite permutation resolution; zero exceedances do not establish zero probability or certainty.

Whole-network structure/global-strength tests and domain-specific tests answer different questions and are reported in separate testing families. Unadjusted omnibus p-values remain exploratory. The 18 domain comparisons use their own BH/Holm correction family. R1 sensitivity and domain comparisons retain their original 1,000 permutations; their non-rejections must be interpreted with the available resolution in mind.

Joint global strength combines symptom-symptom, MRI-MRI and symptom-MRI domains. An omnibus difference alone does not establish a change in symptom-MRI coupling. Individual-edge difference tests do not establish that either separate edge is nonzero or clinically important.

## Sensitivities and equivalence

High correlations between graph weights do not ensure identical significance conclusions. Alternative family-member selections use overlapping samples and are sensitivity analyses rather than independent replications.

Pairwise exploratory TOST results and the network-level claim that all 72 symptom-MRI associations lie within a bound are distinct questions. The all-pairs intersection-union claim requires every pair to pass its unadjusted TOST; separately adjusted pairwise results answer a different multiplicity question. Bounds of ±.05 and ±.03 are exploratory rather than validated clinical thresholds.

Regularised-edge bootstrap intervals are pointwise, conditional on preprocessing, and are not significance tests. Missing-data analyses do not establish MAR or eliminate all selection bias. Network differences do not establish within-person change or causation.
