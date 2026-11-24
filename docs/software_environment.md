# Software environment

The local computational verification was performed with R 4.6.1 (2026-06-24 ucrt), platform x86_64-w64-mingw32/x64, on Windows. The directly checked package versions are in `validated_package_versions.csv`. These versions were observed in the November verification session; this is not a reconstruction of an unobserved historical environment.

`setup_packages.R` installs missing packages only. It does not pin dependencies or force upgrades. A direct-version table is not a complete transitive lockfile. Some packages require system libraries or build tools when suitable binaries are unavailable.

The NCT stage records package versions and session information privately and includes the environment in checkpoint signatures. Changing software, input files or the relevant code after a partial run stops checkpoint reuse. Reconstructed graphs and observed statistics must agree before new p-values are accepted.

For a complete lockfile from an authorised execution environment, run `freeze_environment.R` after verification. It writes `renv.lock`, package versions and session information to the configured external work directory. Review the lockfile before adding it to a versioned code release; do not copy the rest of the work directory. A lockfile describes the captured execution environment, not automatically a historical R1 environment.
