######################################
######################################
#### %>%

#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL


######################################
######################################
#### check...()

#' @title Check that arguments supplied via ... are allowed
#' @description This function checks that arguments supplied via ... are allowed. This function was written to support other functions, specifically via the return of a helpful error message if arguments that cannot be supplied via ... have been supplied. The function is not intended for general use.
#'
#' @param not_allowed A character vector of the names of function arguments that are not allowed.
#' @param ... Other arguments
#'
#' @return The function checks other arguments supplied via ...; if these contain an argument that is not allowed, the function returns an error. Otherwise, nothing is returned.
#'
#' @author Edward Lavender
#' @keywords "internal"

check... <- function(not_allowed,...){
  l <- list(...)
  if(any(names(l) %in% not_allowed)){
    trouble <- names(l)[names(l) %in% not_allowed]
    msg <- paste0("Additional arguments (", paste(trouble, collapse = ", "),
                  ") have been passed to the function via ... which are implemented internally or need to be supplied via other function arguments. Implement these options via appropriate function arguments, if possible, or do not supply them.")
    stop(msg)
  }
}

######################################
######################################
#### check_value()

#' @title Check the input value to a parent function argument
#' @description Within a function, this function checks the value of an input to an argument of that function. If the input value is supported, the function simply returns this value. If the input is not supported, the function either throws an error or returns a warning and the default value.
#'
#' @param arg A character string which defines the argument of the parent function.
#' @param input The input to an argument of a parent function.
#' @param supp A vector of supported input values for the argument in the parent function.
#' @param warn A logical input which defines whether or not to return a warning and the default value (see \code{default}) or an error.
#' @param default The default input value for the parent function.
#'
#' @return The function returns \code{input}, an error or a warning and \code{default} depending on whether or not \code{input} is within \code{supp} (i.e., whether or not the input to the argument of a parent function is supported).
#'
#' @author Edward Lavender
#' @keywords internal
#'

check_value <- function(arg = deparse(substitute(input)), input, supp, warn = TRUE, default = supp[1]){
  # If the input is not in a vector of supported arguments...
  if(!(input %in% supp)){
    ## Provide a warning and revert to the default
    if(is.character(input)) input <- paste0("'", input, "'")
    if(warn){
      if(is.character(default)) default <- paste0("'", default, "'")
      warning(paste0("Argument '", arg, "' = ", input, " is not supported; defaulting to ", arg, " = ", default, ".\n"))
      input <- default
    } else{
      if(is.character(supp)) supp <- paste0("'", supp, "'")
      stop(paste0("Argument '", arg, "' = ", input, " is not supported. Supported option(s): ", paste0(supp, collapse = ", "), "."))
    }
  }
  # Return input
  return(input)
}


###################################
###################################
#### check_class()

#' @title Check the class of an function input to a parent function
#' @description This function checks that the class of an input to a parent function is appropriate. If not, the function either produces a helpful error message or returns a warning.
#' @param arg A character string which defines the argument of the parent function.
#' @param input The input to an argument of a parent function.
#' @param if_class (optional) A character vector of classes of object. If supplied, the function will only proceed to check the class of the object if the \code{class(input)} is one of \code{if_class}. This is useful if \code{check_class()} is implemented in a loop.
#' @param to_class The required class of the input.
#' @param type A character which specifies whether to return an error (\code{"stop"}) or a warning ("warning").
#' @param coerce_input A function used to coerce \code{input} to the correct object type, if \code{type = "warning"}.
#' @return The function checks the class of the input. If the class is not the same as required by the parent function (i.e., as specified by \code{class}), the function returns a helpful error message, or a warning and an object whose class has been coerced to the correct class.
#'
#' @author Edward Lavender
#' @keywords internal
#'

check_class <-
  function(arg = deparse(substitute(input)), input, if_class = NULL, to_class, type = "stop", coerce_input){

    #### Define whether or not to proceed:
    # Only proceed if if_class is NULL or, if supplied, then only proceed if the class of the object
    # ... is of type(s) in if_class
    proceed <- FALSE
    if(is.null(if_class)){
      proceed <- TRUE
    } else{
      if(inherits(input, if_class)) proceed <- TRUE
    }

    #### Check the class, if required
    if(proceed){
      # If the object is not of the necessary class
      if(!inherits(input, to_class)){
        # Either stop...
        if(type == "stop"){
          msg <- paste0("Argument '", arg, "' must be of class '", to_class, "', not class(es): '", paste(class(input), collapse = "', '"), "'.")
          stop(msg)
          # Or print a warning and use coerce_input() to convert the object to the desired class.
        } else if(type == "warning"){
          msg <- paste0("Argument '", arg, "' coerced to class '", to_class, "' from class(es): '", paste(class(input), collapse = "', '"), "'.")
          warning(msg)
          input <- coerce_input(input)
        }
      }
    }

    #### If we've passed all checks, return the input (possibly coerced to a new class)
    return(input)
  }


######################################
######################################
#### check_names()

#' @title Check the names of an object contain required names
#' @description This function checks whether required names are contained within an object. If the object does not contain any/all required names (the precise criteria is controlled by the user), the function returns a helpful error message.
#' @param arg A character string which defines the argument of the parent function.
#' @param input An object for which the names need to be checked.
#' @param req A character vector of required names.
#' @param extract_names A function which is used to extract names from \code{input}, such as \code{\link[base]{names}} or \code{\link[base]{colnames}}.
#' @param type A function which defines the failure criteria. For example, if \code{type = all}, the function will return an error unless all the names in \code{req} are contained within \code{input}. This is the default. If \code{type = any}, the function will return an error only if none of the names in \code{req} are contained within \code{input}.
#' @return If the input fails the check, the function returns a helpful error message. Otherwise, nothing is returned.
#'
#' @author Edward Lavender
#' @keywords internal
#'

check_names <- function(arg = deparse(substitute(input)), input, req, extract_names = names, type = any){
  input_names <- extract_names(input)
  if(!type(req %in% input_names)){
    req_names_missing <- req[which(!(req %in% input_names))]
    msg <- paste0("Argument ", arg, " does not contain ", deparse(substitute(type)), " required names. The following name(s) are missing:",
                  paste0("'", req_names_missing, collapse = ", "),
                  "'.")
    stop(msg)
  }
}


######################################
######################################
#### check_tz()

#' @title Check the timezone of an object and force UTC if absent
#' @description This function checks the time zone of an inputted  object. If the object is of class Date or POSIXct and a time zone is absent, then "UTC" is forced. Otherwise, the object is returned unchanged.
#' @param arg (optional) A character string which defines the argument of the parent function.
#' @param input An object.
#' @return An object as inputted in which, if the object is of class Date or POSIXct and a time zone is absent, time zone "UTC" is forced.
#'
#' @author Edward Lavender
#' @keywords internal

check_tz <-
  function(arg = deparse(substitute(input)), input){
    if(inherits(input, "Date") | inherits(input, "POSIXct")){
      if(lubridate::tz(input) == ""){
        msg <- paste0("Argument '", arg, "' time zone currently ''; tz forced to UTC.")
        warning(msg)
        lubridate::tz(input) <- "UTC"
      }
    }
    return(input)
  }


###################################
###################################
#### check_named_list()

#' @title Check that a list is named
#' @description This function checks that the top level of a list is named (ignoring empty lists if requested). If the list is not named, the function returns a helpful error message. Otherwise, the list is returned unchanged. This is particularly useful within functions that use \code{\link[base]{do.call}} to evaluate lists of arguments.
#' @param arg (optional) A character string which defines the argument of a parent function.
#' @param input A list.
#' @param ignore_empty A logical input which defines whether or not to ignore empty lists.
#' @return The function returns a helpful error message for unnamed lists (ignoring empty lists if requested) or the inputted list unchanged.

#' @author Edward Lavender
#' @keywords internal

check_named_list <- function(arg = deparse(substitute(input)), input, ignore_empty = TRUE){
  if(!any("list" %in% class(input))) stop(paste0("Argument '", arg, "' must be of class list."))
  list_is_empty <- (length(input) == 0)
  if(!list_is_empty | !ignore_empty){
    if(is.null(names(input)) | any(names(input) %in% "")){
      msg <- paste0("Argument '", arg, "' must be a named list.")
      stop(msg)
    }
  }
  return(input)
}


######################################
######################################
#### check_dir()

#' @title Check a directory exists
#' @description This function checks whether a directory exists and, if not, returns an informative error message. The inputted directory can be edited with the addition of a '/' if requested.
#' @param arg (optional) A character string which defines the argument of a parent function.
#' @param input A character string which defines a directory.
#' @param check_slash A logical input that defines whether or not to check the end of a string for '/'. If \code{TRUE} and a '/' is lacking, this is added to the returned directory.
#' @return The function checks whether or not a directory exists. If so, the function returns either the directory as inputted, or the directory with a '/' added to the end. If not, the function returns an informative error message.
#'
#' @author Edward Lavender
#' @keywords internal
#'

check_dir <- function(arg = deparse(substitute(input)),
                      input,
                      check_slash = FALSE){

  #### Check the directory exists
  if(!dir.exists(input)){
    stop(paste0("The directory inputted to the argument '", arg, "' ('", input, "') does not exist."), call. = FALSE)
  }

  #### Check the directory ends in a /
  if(check_slash){
    end_is_slash <- substr(input, nchar(input), nchar(input)) == "/"
    if(!end_is_slash){

      message(paste0("'/' added to the directory inputted to the argument '", arg, "' ('", input, "')."))
      input <- paste0(input, "/")
    }
  }

  #### Return input, possibly updated with / if checks passed
  return(input)
}


#### End of code.
######################################
######################################
