## copied from glmmTMB
## action: message, warning, stop
check_dots <- function(..., action="stop") {
    L <- list(...)
    if (length(L)>0) {
        FUN <- get(action)
        FUN("unknown arguments: ",
            paste(names(L), collapse=","))
    }
    return(NULL)
}
