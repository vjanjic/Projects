@SET_MAKE@

bin_PROGRAMS = bodytrack_seq
bodytrack_SOURCES = DMatrix.h \
                    GeneralUtil.h \
                    ParticleFilter.h \
                    SmallVectors.h \
                    Vector3.h \
                    system.h \
                    AnnealingFactor.h \
                    AnnealingFactor.cpp \
                    BodyGeometry.h \
                    BodyGeometry.cpp \
                    BodyPose.h \
                    BodyPose.cpp\
                    CameraModel.h \
                    CameraModel.cpp \
                    CovarianceMatrix.h \
                    CovarianceMatrix.cpp \
                    ImageMeasurements.h \
                    ImageMeasurements.cpp\
                    ImageProjection.h \
                    ImageProjection.cpp \
                    RandomGenerator.h \
                    RandomGenerator.cpp \
                    TrackingModel.h \
                    TrackingModel.cpp \
                    main.cpp

AM_CPPFLAGS = -I$(VPATH)/../FlexImageLib
bodytrack_LDADD = ../FlexImageLib/libflex.la -lm



