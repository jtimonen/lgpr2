# Define the purely abstract base class Transform
Transform <- R6::R6Class("Transform",
  lock_class = TRUE,
  private = list(
    suffix = NULL
  ),
  public = list(
    print = function() {
      cat("An R6 object of class ", class(self)[1], ".\n", sep = "")
    },
    forward = function(x) {
      stop("'forward()' should be implemented by inheriting class")
    },
    backward = function(x) {
      stop("'backward()' should be implemented by inheriting class")
    },
    add_suffix = function(str) {
      checkmate::assert_character(str)
      if (is.null(private$suffix)) {
        return(str)
      }
      paste0(str, "_", private$suffix)
    }
  )
)

# Define the LinearTransform (shifts and scales data linearly)
LinearTransform <- R6::R6Class("LinearTransform",
  lock_class = TRUE,
  inherit = Transform,
  private = list(
    multiplier = 1.0,
    offset = 0.0,
    suffix = "shifted_and_scaled"
  ),
  public = list(
    initialize = function(multiplier = 1.0, offset = 0.0) {
      private$multiplier <- multiplier
      private$offset <- offset
    },
    forward = function(x) {
      (x + private$offset) * private$multiplier
    },
    backward = function(x) {
      x / private$multiplier - private$offset
    }
  )
)

# Define the IdentityTransform (doesn't do anything)
IdentityTransform <- R6::R6Class("IdentityTransform",
  lock_class = TRUE,
  inherit = LinearTransform,
  private = list(
    suffix = NULL
  ),
  public = list(
    initialize = function() {
      super$initialize(multiplier = 1.0, offset = 0.0)
    }
  )
)


# Define the ScaleTransform (scales data linearly)
ScaleTransform <- R6::R6Class("ScaleTransform",
  lock_class = TRUE,
  inherit = LinearTransform,
  private = list(
    suffix = "scaled"
  ),
  public = list(
    initialize = function(multiplier = 1.0) {
      super$initialize(multiplier = multiplier)
    }
  )
)

# Define the MaxScaleTransform (like ScaleTransform but can be set_using_data)
MaxScaleTransform <- R6::R6Class("MaxScaleTransform",
  lock_class = TRUE,
  inherit = ScaleTransform,
  public = list(
    set_using_data = function(x_data) {
      checkmate::assert_numeric(x_data)
      max_x <- max(x_data)
      checkmate::assert_number(max_x, lower = 1e-12)
      mult <- 1.0 / max_x
      MaxScaleTransform$new(multiplier = mult)
    }
  )
)


# Define the UnitScaleTransform (like ScaleTransform but can be set_using_data)
UnitScaleTransform <- R6::R6Class("UnitScaleTransform",
  lock_class = TRUE,
  inherit = LinearTransform,
  private = list(
    suffix = "unit"
  ),
  public = list(
    # set so that x_data will map to interval [-1, 1]
    set_using_data = function(x_data) {
      checkmate::assert_numeric(x_data)
      x1 <- max(x_data)
      x2 <- min(x_data)
      ofs <- -(x1 + x2) / 2
      mult <- 2 / (x1 - x2)
      UnitScaleTransform$new(multiplier = mult, offset = ofs)
    }
  )
)
