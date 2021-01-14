cmake_minimum_required(VERSION 3.17)
project(klevel)

set(CMAKE_CXX_STANDARD 14)
add_definitions(-DILOUSESTL -DIL_STD)

include_directories(
        header
        utk
        /opt/ibm/ILOG/CPLEX_Studio1210/cplex/include/
        /opt/ibm/ILOG/CPLEX_Studio1210/concert/include/
        /home/jiahao/library/qhull-2020.2/src/libqhull_r/
        /home/jiahao/disk/Projects/klevel/lpsolver/
)

add_executable(klevel
        main.cpp

        /home/jiahao/library/qhull-2020.2/src/libqhull_r/libqhull_r.h
        /home/jiahao/library/qhull-2020.2/src/libqhull_r/geom_r.h
        /home/jiahao/library/qhull-2020.2/src/libqhull_r/io_r.h
        /home/jiahao/library/qhull-2020.2/src/libqhull_r/mem_r.h
        /home/jiahao/library/qhull-2020.2/src/libqhull_r/merge_r.h
        /home/jiahao/library/qhull-2020.2/src/libqhull_r/poly_r.h
        /home/jiahao/library/qhull-2020.2/src/libqhull_r/qhull_ra.h
        /home/jiahao/library/qhull-2020.2/src/libqhull_r/qset_r.h
        /home/jiahao/library/qhull-2020.2/src/libqhull_r/random_r.h
        /home/jiahao/library/qhull-2020.2/src/libqhull_r/stat_r.h
        /home/jiahao/library/qhull-2020.2/src/libqhull_r/user_r.h
        /home/jiahao/library/qhull-2020.2/src/libqhull_r/geom2_r.c
        /home/jiahao/library/qhull-2020.2/src/libqhull_r/geom_r.c
        /home/jiahao/library/qhull-2020.2/src/libqhull_r/global_r.c
        /home/jiahao/library/qhull-2020.2/src/libqhull_r/io_r.c
        /home/jiahao/library/qhull-2020.2/src/libqhull_r/libqhull_r.c
        /home/jiahao/library/qhull-2020.2/src/libqhull_r/mem_r.c
        /home/jiahao/library/qhull-2020.2/src/libqhull_r/merge_r.c
        /home/jiahao/library/qhull-2020.2/src/libqhull_r/poly2_r.c
        /home/jiahao/library/qhull-2020.2/src/libqhull_r/poly_r.c
        /home/jiahao/library/qhull-2020.2/src/libqhull_r/qset_r.c
        /home/jiahao/library/qhull-2020.2/src/libqhull_r/random_r.c
        /home/jiahao/library/qhull-2020.2/src/libqhull_r/rboxlib_r.c
        /home/jiahao/library/qhull-2020.2/src/libqhull_r/stat_r.c
        /home/jiahao/library/qhull-2020.2/src/libqhull_r/usermem_r.c
        /home/jiahao/library/qhull-2020.2/src/libqhull_r/userprintf_r.c
        /home/jiahao/library/qhull-2020.2/src/libqhull_r/userprintf_rbox_r.c
        /home/jiahao/library/qhull-2020.2/src/libqhull_r/user_r.c


        header/header.h
        header/klevel.h
        klevel.cpp
        header/hyperplane.h
        hyperplane.cpp

        utk/point.h
        utk/point.cpp
        utk/hypercube.h
        utk/hypercube.cpp
        utk/rnode.h
        utk/rnode.cpp
        utk/rentry.h
        utk/rentry.cpp
        utk/rtree.h
        utk/rtree.cpp
        utk/collection.h
        utk/collection.cpp
        utk/filemem.h
        utk/filemem.cpp
        utk/tgs.h
        utk/tgs.cpp
        utk/object.h
        utk/object.cpp
        utk/skyline.h
        utk/skyline.cpp
        utk/virtualRNode.h
        utk/virtualRNode.cpp
        utk/global.h
        utk/global.cpp
        utk/cellTree.h
        utk/cellTree.cpp
        utk/utk.h
        utk/utk.cpp

        lpsolver/lp_Hash.h
        lpsolver/lp_lib.h
        lpsolver/lp_matrix.h
        lpsolver/lp_mipbb.h
        lpsolver/lp_SOS.h
        lpsolver/lp_utils.h
        )

find_library(ILOCPLEX ilocplex /opt/ibm/ILOG/CPLEX_Studio1210/cplex/lib/x86-64_linux/static_pic)
find_library(CPLEX cplex /opt/ibm/ILOG/CPLEX_Studio1210/cplex/lib/x86-64_linux/static_pic)
find_library(CONCERT concert /opt/ibm/ILOG/CPLEX_Studio1210//concert/lib/x86-64_linux/static_pic/)

find_library(LPSOLVER liblpsolve55.a /home/jiahao/disk/Projects/klevel/lpsolver/)

target_link_libraries(klevel PUBLIC
        ${ILOCPLEX}
        ${CONCERT}
        ${CPLEX}
        libm.so # should use dynamic lib
        libpthread.so # should use dynamic lib
        /home/jiahao/disk/Projects/klevel/lpsolver/liblpsolve55.so
        ${LPSOLVER}
        ${CMAKE_DL_LIBS}
        )

