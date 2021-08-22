#!/usr/bin/bash

t_flag=false
n_flag=false
s_flag=false

regex="^[0-9]+$"

print_usage() {
  printf "Usage: $0 [-t #num] [-s <secuencia>] [-n #num]\n"
  exit 1
}

is_number() {
	if ! [[ "$1" =~ $regex ]]
	then
		printf "Error: -$2 argument must be a number\n"
		print_usage
		exit 1
	fi
}

verify_flags() {
	for flag in "$@"
	do
		if [ $flag = false ]
		then
			printf "Missing mandatory arguments\n"
			print_usage
			exit 1
		fi
	done
}

# Flag configuration
while getopts t:s:n: flag
do
	case $flag in
		t) t_flag=true
		   is_number "$OPTARG" "$flag"
		   if [ $OPTARG \> 0 ]
		   then
		   	TRIES=$OPTARG
		   else
		   	printf "Argument for -t must be greathen than zero\n"
		   	exit 1
		   fi ;;
		   
		s) s_flag=true
		   if [ $OPTARG = random ] || [ $OPTARG = inv ] || [ $OPTARG = sorted ]
		   then
		   	SEQUENCE=$OPTARG
		   else
		   	printf "Invalid identifier: $OPTARG\n"
		   	printf "Identifiers are: random, inv and sorted\n"
		   	print_usage
		   fi ;;

		n) n_flag=true
		   is_number "$OPTARG" "$flag"
		   NUMBER=$OPTARG ;;

		[?]) print_usage ;;

	esac
done

verify_flags $t_flag $n_flag $s_flag

java -jar TestSort.jar $TRIES $SEQUENCE $NUMBER