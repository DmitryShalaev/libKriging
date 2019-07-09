#!/usr/bin/env bash
set -eo pipefail

if [[ "$DEBUG_CI" == true ]]; then
  set -x
fi

# Default configuration when used out of travis-ci
MODE=${MODE:-Debug}
EXTRA_CMAKE_OPTIONS=${EXTRA_CMAKE_OPTIONS:-}
DO_TEST=${DO_TEST:-true}

if [[ -n ${TRAVIS_BUILD_DIR:+x} ]]; then
    cd ${TRAVIS_BUILD_DIR}
fi

mkdir -p build
cd build
cmake \
  -G "Unix Makefiles" \
  -DCMAKE_BUILD_TYPE=${MODE} \
  ${EXTRA_CMAKE_OPTIONS} \
  ..

if [[ "$DO_TEST" == "true" ]]; then
    cmake --build .

    if [[ "$MODE" == "Coverage" ]]; then
        cmake --build . --target coverage --config Coverage
    else
        ctest -C ${MODE} # --verbose
    fi

    cmake --build . --target install --config ${MODE}
    # Show installation directory content
    # find installed
else
    # faster install target if tests are not required
    cmake --build . --target install.lib
fi




