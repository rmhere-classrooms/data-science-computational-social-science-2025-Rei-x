library(shiny)
library(bslib)
library(igraph)

seed_method_options <- c(
  "Losowy wybór 5% węzłów" = "random",
  "Top 5% według stopnia wychodzącego" = "top_out_degree",
  "Top 5% według betweenness" = "top_betweenness",
  "Top 5% według closeness" = "top_closeness",
  "Top 5% według ważonego zasięgu dwuskokowego" = "top_two_hop"
)

sanitize_scores <- function(scores) {
  scores[!is.finite(scores)] <- 0
  scores
}

compute_two_hop_reach <- function(graph, out_deg, edge_weights) {
  ends_mat <- ends(graph, es = E(graph), names = FALSE)
  scores <- numeric(vcount(graph))

  for (idx in seq_len(nrow(ends_mat))) {
    from <- ends_mat[idx, 1]
    to <- ends_mat[idx, 2]
    w <- edge_weights[idx]
    scores[from] <- scores[from] + w * out_deg[to]
  }

  scores
}

compute_seed_scores <- function(graph) {
  edge_weights <- E(graph)$weight
  edge_weights[is.na(edge_weights)] <- 0
  distance_weights <- ifelse(edge_weights > 0, 1 / pmax(edge_weights, 1e-6), Inf)

  out_deg <- degree(graph, mode = "out")

  list(
    top_out_degree = out_deg,
    top_betweenness = sanitize_scores(betweenness(graph, directed = TRUE, weights = distance_weights)),
    top_closeness = sanitize_scores(closeness(graph, mode = "out", weights = distance_weights, normalized = TRUE)),
    # W2HR (Weighted Two-Hop Reach) premiuje węzły, które łączą wysoką wagę wyjścia z sąsiadami posiadającymi duży stopień wychodzący,
    # przez co z większym prawdopodobieństwem uruchamiają szeroką kaskadę w dwóch krokach.
    top_two_hop = sanitize_scores(compute_two_hop_reach(graph, out_deg, edge_weights))
  )
}

select_seed_nodes <- function(graph, method, fraction = 0.05, score_cache = NULL) {
  n <- vcount(graph)
  seed_count <- max(1, ceiling(fraction * n))
  vertices <- seq_len(n)

  if (method == "random") {
    sample(vertices, seed_count)
  } else {
    if (is.null(score_cache) || is.null(score_cache[[method]])) {
      score_cache <- compute_seed_scores(graph)
    }

    if (is.null(score_cache[[method]])) {
      stop("Brak zdefiniowanego sposobu wyboru dla metody: ", method)
    }

    scores <- score_cache[[method]]
    scores[is.na(scores)] <- -Inf
    ordered <- order(scores, decreasing = TRUE, na.last = NA)
    ordered[seq_len(seed_count)]
  }
}

run_independent_cascade_once <- function(graph, seeds, prob_multiplier = 1.0) {
  n <- vcount(graph)
  activated <- rep(FALSE, n)
  seeds <- unique(seeds)
  activated[seeds] <- TRUE

  frontier <- seeds
  counts <- c(sum(activated))
  attempted <- vector("list", n)

  while (length(frontier) > 0) {
    next_frontier <- integer(0)

    for (node in frontier) {
      out_edges <- incident(graph, v = node, mode = "out")
      if (length(out_edges) == 0) {
        next
      }

      edge_pairs <- ends(graph, es = out_edges, names = FALSE)
      targets <- as.integer(edge_pairs[, 2])
      probs <- E(graph)[out_edges]$weight
      probs[is.na(probs)] <- 0
      probs <- probs * prob_multiplier
      probs <- pmin(pmax(probs, 0), 1)
      attempted_neighbors <- attempted[[node]]

      for (idx in seq_along(targets)) {
        neighbor <- targets[idx]

        if (!is.null(attempted_neighbors) && neighbor %in% attempted_neighbors) {
          next
        }

        attempted_neighbors <- c(attempted_neighbors, neighbor)

        if (!activated[neighbor] && runif(1) <= probs[idx]) {
          activated[neighbor] <- TRUE
          next_frontier <- c(next_frontier, neighbor)
        }
      }

      attempted[[node]] <- attempted_neighbors
    }

    if (length(next_frontier) == 0) {
      break
    }

    frontier <- unique(next_frontier)
    counts <- c(counts, sum(activated))
  }

  counts
}

average_cascade_curve <- function(curves) {
  max_len <- max(vapply(curves, length, integer(1)))
  mat <- matrix(NA_real_, nrow = length(curves), ncol = max_len)

  for (i in seq_along(curves)) {
    mat[i, seq_along(curves[[i]])] <- curves[[i]]
  }

  colMeans(mat, na.rm = TRUE)
}

simulate_cascade_method <- function(graph, method, trials, score_cache = NULL, prob_multiplier = 1.0) {
  stopifnot(trials > 0)
  trial_results <- vector("list", trials)

  for (i in seq_len(trials)) {
    seeds <- select_seed_nodes(graph, method, score_cache = score_cache)
    trial_results[[i]] <- run_independent_cascade_once(graph, seeds, prob_multiplier)
  }

  avg_curve <- average_cascade_curve(trial_results)
  new_activations <- c(avg_curve[1], diff(avg_curve))

  data.frame(
    iteration = seq_along(avg_curve) - 1,
    avg_activated = avg_curve,
    new_activations = new_activations,
    method_key = method,
    stringsAsFactors = FALSE
  )
}

ui <- page_sidebar(
  title = "Graf sieci emailowej",
  sidebar = sidebar(
    h4("Informacje o grafie"),
    textOutput(outputId = "graphStats"),
    hr(),
    h4("Symulacja dyfuzji"),
    checkboxGroupInput(
      inputId = "seedMethods",
      label = "Sposób wyboru węzłów startowych",
      choices = seed_method_options,
      selected = seed_method_options
    ),
    sliderInput(
      inputId = "probMultiplier",
      label = "Mnożnik prawdopodobieństwa wij (%)",
      min = 10,
      max = 200,
      value = 100,
      step = 10,
      post = "%"
    ),
    # Iteracje = liczba powtórzeń eksperymentu (1-50) zgodnie z pkt 10 PDF.
    sliderInput(
      inputId = "trialCount",
      label = "Liczba iteracji (powtórzeń)",
      min = 1,
      max = 50,
      value = 10,
      step = 1
    ),
    actionButton("runSim", "Uruchom symulację", class = "btn-primary")
  ),
  plotOutput(outputId = "networkPlot"),
  plotOutput(outputId = "cascadePlot"),
  tableOutput(outputId = "cascadeSummary")
)



server <- function(input, output) {
  # Reactive: wczytanie danych, policzenie wag (wij = cntij / cnti), budowa grafu
  email_graph_reactive <- reactive({
    # Wczytaj oryginalny DF krawędzi (from, to)
    email_data <- read.table(
      "https://bergplace.org/share/out.radoslaw_email_email",
      skip = 2,
      colClasses = c("integer", "integer", "NULL", "NULL")
    )
    colnames(email_data) <- c("from", "to")

    # Zlicz cntij dla par (from, to)
    counts_df <- aggregate(
      x = list(cntij = rep(1L, nrow(email_data))),
      by = list(from = email_data$from, to = email_data$to),
      FUN = sum
    )

    totals_df <- aggregate(cntij ~ from, data = counts_df, sum)
    names(totals_df)[2] <- "cnti"

    counts_df <- merge(counts_df, totals_df, by = "from", all.x = TRUE)
    counts_df$weight <- counts_df$cntij / counts_df$cnti

    g <- graph_from_data_frame(counts_df[, c("from", "to", "weight")], directed = TRUE)
    g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
    stopifnot(vcount(g) == 167, ecount(g) == 5783)

    g
  })

  output$networkPlot <- renderPlot({
    email_graph <- email_graph_reactive()

    edge_w <- 1 + 6 * E(email_graph)$weight

    plot(
      email_graph,
      vertex.size = 5,
      vertex.label = NA,
      edge.width = edge_w,
      edge.arrow.size = 0.3,
      edge.color = "#007bc2",
      vertex.color = "#ff6b6b",
      main = "Graf skierowany sieci emailowej (krawędzie ważone)"
    )
  })

  simulation_results <- eventReactive(input$runSim, {
    methods_selected <- input$seedMethods
    req(length(methods_selected) > 0)

    email_graph <- email_graph_reactive()
    trials <- as.integer(input$trialCount)
    prob_multiplier <- input$probMultiplier / 100
    score_cache <- compute_seed_scores(email_graph)

    method_results <- lapply(methods_selected, function(method_key) {
      result_df <- simulate_cascade_method(email_graph, method_key, trials, score_cache, prob_multiplier)
      label <- names(seed_method_options)[match(method_key, seed_method_options)]
      result_df$method_label <- label
      result_df
    })

    do.call(rbind, method_results)
  }, ignoreNULL = FALSE)

  output$cascadePlot <- renderPlot({
    res <- simulation_results()
    req(!is.null(res), nrow(res) > 0)

    max_iter <- max(res$iteration)
    max_new <- max(res$new_activations, na.rm = TRUE)
    methods <- unique(res$method_label)
    colors <- c("#ff6b6b", "#007bc2", "#1dd1a1", "#ffa502", "#2e86de")

    plot(
      NA,
      xlim = c(0, max_iter),
      ylim = c(0, max_new),
      xlab = "Nr iteracji",
      ylab = "Liczba nowo aktywowanych węzłów",
      main = "Proces dyfuzji informacji - nowe aktywacje na iterację",
      xaxt = "n"
    )
    axis(1, at = pretty(seq(0, max_iter, length.out = min(max_iter + 1, 10))))

    for (i in seq_along(methods)) {
      method_df <- res[res$method_label == methods[i], ]
      lines(method_df$iteration, method_df$new_activations, col = colors[i], lwd = 2)
    }

    legend("topright", legend = methods, col = colors[seq_along(methods)], lwd = 2, bty = "n")
  })

  output$cascadeSummary <- renderTable({
    res <- simulation_results()
    req(!is.null(res), nrow(res) > 0)
    total_nodes <- vcount(email_graph_reactive())

    summary_df <- do.call(rbind, lapply(split(res, res$method_label), function(df) {
      final_row <- df[nrow(df), ]
      data.frame(
        Metoda = df$method_label[1],
        `Śr. aktywne (koniec)` = final_row$avg_activated,
        `Śr. % aktywnych` = 100 * final_row$avg_activated / total_nodes,
        `Liczba iteracji` = nrow(df) - 1,
        check.names = FALSE
      )
    }))

    summary_df
  }, digits = 1, striped = TRUE, bordered = TRUE)

  output$graphStats <- renderText({
    email_graph <- email_graph_reactive()

    v <- vcount(email_graph)
    e <- ecount(email_graph)

    paste0(
      "Wierzchołki: ", v, "\n",
      "Krawędzie: ", e
    )
  })
}

shinyApp(ui = ui, server = server)
