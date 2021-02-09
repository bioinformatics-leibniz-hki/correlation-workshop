#!/usr/bin/env R

#
# Initialization script called by Rstudio Server
#

# Auto start project
setHook("rstudio.sessionInit", function(newSession) {
  if (newSession && is.null(rstudioapi::getActiveProject())) {
    rstudioapi::openProject(".")
  }
}, action = "append")
