#### Package Startup ####

.onAttach <- function(
    libname,
    pkgname
) {

    version <- utils::packageVersion(
        pkgname
    )

    packageStartupMessage(
        paste0(
            "EconR ", version, "\n",
            "Econometric models with consistent estimation semantics\n",
            "David Bustamante Lazo\n",
            "GitHub: https://github.com/DavidBustamanteL/EconR"
        )
    )
}