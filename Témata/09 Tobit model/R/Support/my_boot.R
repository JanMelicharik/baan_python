# Vazeny bootstrap
# index ... S*1 vektor bootstrapovych indexu
# y ... N x 1 vektor bootstrapovaneho vzorku
# S ... velikost bootstrapoveho vyberu
# w ... bootstrapove vahy (optional)

my_boot <- function(y,S,w = NULL) {
    N <- length(y)

    if(is.null(w)) {
        w <- 1/N
    }

    # vektor indexu
    boot_index <- mat.or.vec(S,1)

    # vytvoreni bodu intervalu, kde sirka dilcich intervalu reprezentuje
    # stanovenou vahu
    w_int <- cumsum(w)

    for(s in 1:S) {
        # runif*w_int(end) ...nahodne cislo z uniformniho rozdeleni na intervalu
        # (0;w_int(length(w)) a ziskani indexu jako suma intervalu splnujici
        # podminku, ze hodnota krajni meze intervalu je mensi nez toto nahodne
        # vygenerovane cislo
        boot_index[s] <- sum(w_int < (runif(1)*w_int[length(w)])) + 1
    }
    return(boot_index)
}
