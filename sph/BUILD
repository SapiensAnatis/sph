# Library excludes main.cpp and is exposed to unit testing

cc_library(
    name = "sph-lib",
    srcs = glob(["*.cpp"]),
    hdrs = glob(["*.hpp"]),
    visibility = ["//unittest:__pkg__"],
)

cc_binary(
    name = "sph",
    srcs = glob(["*.cpp"]),
    deps = [
        ":sph-lib",
    ]
)