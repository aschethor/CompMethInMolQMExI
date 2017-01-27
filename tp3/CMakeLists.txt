set(SOURCE_FILES
        md.c
        force.c
        init.c
        integrate.c
        mdloop.c
        ran_uniform.c
        ran_uniform.h
        readdat.c
        sample_diff.c
        sample_rdf.c
        system.c
        system.h
        writepdb.c)
add_executable(tp3 ${SOURCE_FILES})