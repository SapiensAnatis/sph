# Library excludes main.cpp and is exposed to unit testing

cc_library(
    name = "sph-lib",
    srcs = ["setup.cpp"],
    hdrs = ["setup.hpp", "particle.hpp"],
    visibility = ["//unittest:__pkg__"],
)

cc_binary(
    name = "sph",
    srcs = ["*.cpp"],
    deps = [
        ":sph-lib",
    ]
)