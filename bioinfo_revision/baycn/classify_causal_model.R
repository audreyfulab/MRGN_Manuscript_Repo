#' Classify a probabilistic adjacency matrix into one of the causal models
#'
#' Uses only the first three rows/columns of \code{post.adj}, corresponding
#' (in order) to V, T1, T2. Rows are parents, columns are children, so
#' \code{post.adj[i, j]} is Pr(i -> j).
#'
#' Edge existence and direction rules for a pair (i, j):
#'   p_ij  = Pr(i -> j)
#'   p_ji  = Pr(j -> i)
#'   p0    = 1 - p_ij - p_ji                    (Pr(no edge))
#'   present   : p0 < edge.presence
#'   directed  : |p_ij - p_ji| > edge.direction  (direction = larger prob)
#'   undirected: present but not directed
#'
#' @param post.adj A square probabilistic adjacency matrix (or data.frame)
#'   whose first 3 rows/columns correspond to V, T1, T2 in that order.
#' @param edge.presence Threshold on p0 for an edge to be considered present.
#'   Default 0.5.
#' @param edge.direction Threshold on |p_ij - p_ji| for an edge to be
#'   considered directed. Default 0.2.
#'
#' @return A list with:
#'   \item{model}{Character label: "M0", "M1.1", "M1.2", "M2.1", "M2.2",
#'     "M3", "M4", or "Other".}
#'   \item{edges}{Named character vector giving the resolved status of each
#'     of the three pairs (V_T1, V_T2, T1_T2): "none", "A->B", "B->A", or
#'     "A--B" (undirected).}
#'   \item{probs}{Named numeric vector of the six directed-edge posterior
#'     probabilities: V_to_T1, V_from_T1, V_to_T2, V_from_T2, T1_to_T2,
#'     T1_from_T2 (i.e. Pr(V->T1), Pr(V<-T1), Pr(V->T2), Pr(V<-T2),
#'     Pr(T1->T2), Pr(T1<-T2)).}
#'
#' @examples
#' # V->T1->T2 (M1.1)
#' m <- matrix(0, 3, 3, dimnames = list(c("V","T1","T2"), c("V","T1","T2")))
#' m["V","T1"] <- 0.9
#' m["T1","T2"] <- 0.85
#' classify_causal_model(m)
classify_causal_model <- function(post.adj, edge.presence = 0.5, edge.direction = 0.2) {

  nodes <- c("V", "T1", "T2")
  m <- as.matrix(post.adj[1:3, 1:3])
  rownames(m) <- nodes
  colnames(m) <- nodes

  # Resolve the status of a single pair of nodes (a = row-ish, b = col-ish;
  # order only affects the labeling of the returned string, not the logic).
  get_edge <- function(a, b) {
    p_ab <- m[a, b]   # Pr(a -> b)
    p_ba <- m[b, a]   # Pr(b -> a)
    p0   <- 1 - p_ab - p_ba

    if (p0 >= edge.presence) {
      return("none")
    }

    if (abs(p_ab - p_ba) > edge.direction) {
      if (p_ab > p_ba) {
        paste0(a, "->", b)
      } else {
        paste0(b, "->", a)
      }
    } else {
      paste0(a, "--", b)  # present but undirected
    }
  }

  e_V_T1  <- get_edge("V", "T1")
  e_V_T2  <- get_edge("V", "T2")
  e_T1_T2 <- get_edge("T1", "T2")

  probs <- c(
    V_to_T1    = m["V", "T1"],
    V_from_T1  = m["T1", "V"],
    V_to_T2    = m["V", "T2"],
    V_from_T2  = m["T2", "V"],
    T1_to_T2   = m["T1", "T2"],
    T1_from_T2 = m["T2", "T1"]
  )

  model <- "Other"

  if (e_V_T1 == "V->T1" && e_V_T2 == "none" && e_T1_T2 == "none") {
    model <- "M0"
  } else if (e_V_T1 == "none" && e_V_T2 == "V->T2" && e_T1_T2 == "none") {
    model <- "M0"
  } else if (e_V_T1 == "V->T1" && e_V_T2 == "none" && e_T1_T2 == "T1->T2") {
    model <- "M1.1"
  } else if (e_V_T1 == "none" && e_V_T2 == "V->T2" && e_T1_T2 == "T2->T1") {
    model <- "M1.2"
  } else if (e_V_T1 == "V->T1" && e_V_T2 == "none" && e_T1_T2 == "T2->T1") {
    model <- "M2.1"
  } else if (e_V_T1 == "none" && e_V_T2 == "V->T2" && e_T1_T2 == "T1->T2") {
    model <- "M2.2"
  } else if (e_V_T1 == "V->T1" && e_V_T2 == "V->T2" && e_T1_T2 == "none") {
    model <- "M3"
  } else if (e_V_T1 == "V->T1" && e_V_T2 == "V->T2" && e_T1_T2 != "none") {
    model <- "M4"
  }

  list(
    model = model,
    edges = c(V_T1 = e_V_T1, V_T2 = e_V_T2, T1_T2 = e_T1_T2),
    probs = probs
  )
}

