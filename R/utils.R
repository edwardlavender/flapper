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
#### check functions
# Source: utils.add: https://github.com/edwardlavender/utils.add
# 19/06/2020

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
#' @examples
#' #### Example (1) Imagine we have a function in wich xlim and ylim cannot be supplied via ...
#' # Internally, within that function, we can implement check as follows:
#' pf <- function(...){
#'       check...(not_allowed = c("xlim", "ylim"),...)
#'       plot(1:10, 1:10, xlim = c(1, 10), ylim = c(1, 10),...)
#'       }
#' # This works:
#' pf(col = "red")
#' # This returns an error
#' \dontrun{
#' pf(col = "red", xlim = c(1, 15))
#' }
#'
#' @author Edward Lavender
#' @export


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
#' @description Within a function, this function checks the value of an input to an argument of that function. If the input value is supported, the function simply returns this value. If the input is not supported, the function returns a warning and the default value. This function is designed to be implemented internally within functions and not intended for general use.
#'
#' @param arg A character string which defines the argument of the parent function.
#' @param input The input to an argument of a parent function.
#' @param supp A vector of supported input values for the argument in the parent function.
#' @param default The default input value for the parent function.
#'
#' @return The function returns \code{input} or \code{default} (the latter with a warning) depending on whether or not \code{input} is within \code{supp} (i.e., whether or not the input to the argument of a parent function is supported).
#'
#' @examples
#'
#' #### Define an example function:
#' # The function returns 1 or 2, depending on the input to 'output'
#' return_1_or_2 <- function(output = 1){
#'   # Check the output, changing the output to the default if necessary
#'   output <- check_value(arg = "output", input = output, supp = 1:2, default = 1)
#'   # Return a value according to 'output'
#'   if(output == 1) return(1) else if(output == 2) return(2)
#' }
#'
#' #### Example (1): If a supported input to output is provided, everything works perfectly:
#' return_1_or_2(1); return_1_or_2(2)
#'
#' #### Example (2): # If an unsupported input to output is provided,
#' # ... the default output is used with a warning:
#' \dontrun{
#' return_1_or_2(3)
#' }
#'
#' @author Edward Lavender
#' @export
#'

check_value <- function(arg = deparse(substitute(input)), input, supp, default = supp[1]){
  # If the input is not in a vector of supported arguments...
  if(!(input %in% supp)){
    # Provide a warning and revert to the default
    warning(paste0("Argument '", arg, "' = ", input, " is not supported; defaulting to ", arg, " = ", default, ".\n"))
    input <- default
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
#' @examples
#' #### Example (1): Implementation using default options outside of a function
#' # Imagine we have an argument, x, to a function, the input to which must be a list.
#' # This input passes the check:
#' check_class(arg = "x", input = list(), to_class = "list")
#' \dontrun{
#' # This input fails the check:
#' check_class(arg = "x", input = list(), to_class = "Date")
#' }
#'
#' #### Example (2): Implementation within a parent function
#' nest_list_in_list <- function(x){
#'   check_class(arg = "x", input = x, to_class = "list")
#'   if(inherits(x, "list")){
#'     return(list(x))
#'   }
#' }
#' nest_list_in_list(list())
#' \dontrun{
#' nest_list_in_list("a")
#' }
#'
#' #### Example (3) Return a warning rather than an error
#' x <- as.POSIXct("2016-01-01")
#' check_class(arg = "x", input = x, to_class = "Date",
#'                   type = "warning", coerce_input = as.Date)
#'
#' #### Example (4) Only act on objects of a certain class; otherwise, return objects unchanged.
#' # In this case the function checks x:
#' check_class(arg = "x", input = x,
#'                   if_class = c("POSIXct", "Date"),
#'                   to_class = "Date", type = "warning", coerce_input = as.Date)
#' # In this case the function does not check x
#' check_class(arg = "x", input = 5,
#'                   if_class = c("POSIXct", "Date"),
#'                   to_class = "Date", type = "warning", coerce_input = as.Date)
#'
#' @author Edward Lavender
#' @export
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
#' @return If the input fails the check, the function returns a helpful error mesaage. Otherwise, nothing is returned.
#' @examples
#' \dontrun{
#' check_names(input = list(x = 1, y = 1),
#'             req = c("a", "b", "c"),
#'             extract_names = names,
#'             type = all)
#' check_names(input = data.frame(b = 1, a = 1),
#'             req = c("x", "y"),
#'             extract_names = colnames,
#'             type = any)
#' }
#'
#' @author Edward Lavender
#' @export
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
#' @examples
#' check_tz(input = as.POSIXct("2016-01-01"))
#' check_tz(arg = "t", input = as.POSIXct("2016-01-01"))
#' check_tz(arg = "t", input = as.POSIXct("2016-01-01", tz = "UTC"))
#'
#' @author Edward Lavender
#' @export

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
#'
#' @examples
#' # This returns input unchanged:
#' check_named_list(input = list(), ignore_empty = TRUE)
#' check_named_list(input = list(a = "b"))
#' # This returns an error:
#' \dontrun{
#' #' check_named_list(input = list(), ignore_empty = FALSE)
#' }
#' # This returns an error which includes the argument name:
#' \dontrun{
#' check_named_list(arg = "x", input = list(1))
#' }
#'
#' @author Edward Lavender
#' @export

check_named_list <- function(arg = deparse(substitute(input)), input, ignore_empty = TRUE){
  if(plotrix::listDepth(input) > 1){
    warning("Input list of check_named_list() is of depth > 1; only the top level is checked.")
  }
  list_is_empty <- (length(input) == 0)
  if(!list_is_empty | !ignore_empty){
    if(is.null(names(input)) | any(names(input) %in% "")){
      msg <- paste0("Argument '", arg, "' must be a named list.")
      stop(msg)
    }
  }
  return(input)
}


#### End of code.
######################################
######################################

