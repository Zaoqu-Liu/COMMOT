#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2025
#  Saez Lab, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Alberto Valdeolivas
#                  Dénes Türei (turei.denes@gmail.com)
#                  Attila Gábor
#
#  Distributed under the MIT (Expat) License.
#  See accompanying file `LICENSE` or find a copy at
#      https://directory.fsf.org/wiki/License:Expat
#
#  Website: https://r.omnipathdb.org/
#  Git repo: https://github.com/saezlab/OmnipathR
#


#' Keep the curl handle in httr2 under the `handle` attribute of `request`
#'
#' @importFrom utils packageVersion
#' @importFrom logger log_trace
#' @noRd
patch_httr2_keep_handle <- function() {

    PATCHED <- 'omnipathr_patched'

    if (requireNamespace('httr2', quietly = TRUE)) {

        ns <- asNamespace('httr2')
        
        # Check httr2 version for compatibility
        httr2_version <- tryCatch(
            utils::packageVersion('httr2'),
            error = function(e) '0.0.0'
        )
        
        # For httr2 >= 1.0.0, the API changed and patch may not be needed
        if (httr2_version >= '1.0.0') {
            logger::log_trace(
                'httr2 version %s detected, checking if patch is needed.',
                httr2_version
            )
            # Check if req_perform1 exists
            if (!'req_perform1' %in% names(ns)) {
                logger::log_trace('req_perform1 not found in httr2, skipping patch.')
                return(invisible(NULL))
            }
        }
        
        # Try to get req_perform1, it may not exist in newer versions
        original <- tryCatch(
            get('req_perform1', ns),
            error = function(e) {
                logger::log_trace('req_perform1 not found in httr2 namespace.')
                return(NULL)
            }
        )
        
        if (is.null(original)) {
            logger::log_trace('Skipping httr2 patch: req_perform1 not available.')
            return(invisible(NULL))
        }

        if (is.null(attr(original, PATCHED))) {

            # Get the formal arguments of the original function
            original_args <- names(formals(original))
            logger::log_trace('req_perform1 has arguments: %s', 
                            paste(original_args, collapse = ', '))

            # Create patched function that matches the signature
            if ('req_prep' %in% original_args && 'resend_count' %in% original_args) {
                # New httr2 API (>= 1.0.0)
                logger::log_trace('Using new httr2 API signature for patch.')
                patched <- function(req, req_prep, path = NULL, handle = NULL, resend_count = 0) {
                    if (!is.null(handle)) {
                        attr(req, 'handle') <- handle
                    }
                    original(req, req_prep, path = path, handle = handle, resend_count = resend_count)
                }
            } else {
                # Old httr2 API (< 1.0.0)
                logger::log_trace('Using old httr2 API signature for patch.')
                patched <- function(req, path = NULL, handle = NULL) {
                    if (!is.null(handle)) {
                        attr(req, 'handle') <- handle
                    }
                    original(req, path = path, handle = handle)
                }
            }

            attr(patched, PATCHED) <- TRUE

            patch_ns('req_perform1', patched, ns)

        } else {

            patched <- original

        }

        patch_ns('req_perform1', patched, asNamespace('OmnipathR'), add = TRUE)

    }

}
