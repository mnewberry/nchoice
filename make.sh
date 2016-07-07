#!/bin/sh

ocamlbuild -pkgs gsl,pcre,batteries assort.native
ocamlbuild -pkgs gsl,pcre,batteries endpoints.native
ocamlbuild -pkgs gsl,pcre,batteries fixation_mcmc_exact.native
ocamlbuild -pkgs gsl,pcre,batteries fixation_mcmc.native
ocamlbuild -pkgs gsl,pcre,batteries pi.native
