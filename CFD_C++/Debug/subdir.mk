################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../gauss_elimination.cpp \
../gauss_seidel.cpp \
../gs.cpp \
../test.cpp \
../thomas.cpp 

OBJS += \
./gauss_elimination.o \
./gauss_seidel.o \
./gs.o \
./test.o \
./thomas.o 

CPP_DEPS += \
./gauss_elimination.d \
./gauss_seidel.d \
./gs.d \
./test.d \
./thomas.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I/usr/include/python2.7 -I/usr/local/lib64/python2.7/site-packages/numpy/core/include -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


