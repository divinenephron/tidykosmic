# Create a temporary environment to load the
# haemloglobin data into to avoid overwriting
# the existing "haemoglobin" variable.
hemoglobin_env <- new.env()
load(system.file("data", "haemoglobin.rda",
                 package = "tidykosmic"),
     envir = hemoglobin_env)
hemoglobin <- hemoglobin_env$haemoglobin
rm(hemoglobin_env)