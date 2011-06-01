################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../ils.c \
../mils.c \
../qrmcp.c \
../reduction.c \
../search.c \
../test.c 

OBJS += \
./ils.o \
./mils.o \
./qrmcp.o \
./reduction.o \
./search.o \
./test.o 

C_DEPS += \
./ils.d \
./mils.d \
./qrmcp.d \
./reduction.d \
./search.d \
./test.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


