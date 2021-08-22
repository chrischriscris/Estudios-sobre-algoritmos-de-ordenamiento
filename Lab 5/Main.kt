import kotlin.system.measureTimeMillis
import kotlin.math.*
import kotlin.system.exitProcess
import kotlin.Double.Companion.POSITIVE_INFINITY
import kotlin.Double.Companion.NEGATIVE_INFINITY

/**
 * Prints an error [message], indicates the usage of the
 * command line and halts the execution of the program.
 */
fun error(message: String) {
	println(message)
	println("Usage: ./runSorlib.sh [-t #num] [-s <secuencia>] [-n #n_1 ... #n_m] [-a <alg>] [-g <figura>]")
	exitProcess(1)
}

/**
 * Parse the [args] array and extracts the data which
 * will be used in the main() function regardless of
 * the order of the flags in the array.
 *
 * @return A triple with the number of tries,
 * an array with the type of sequence, the set
 * of algorithms and the filename, and an array
 * of the sizes of the sequence.
 */
fun parseCommands(args: Array<String>):
    Triple<Int, Array<String>, Sequence> {
	val validFlags = arrayOf("-t", "-s", "-n", "-a", "-g")

	// Return variables
	var tries = 0
	val parsedArgs = Array<String>(3, { "" })
	val seqSizes: Array<Int?> = arrayOfNulls(args.size)

	// Check every mandatory argument was provided
	for ((i, flag) in validFlags.withIndex()) {
		if (i != 4 && flag !in args) {
			error("Missing mandatory arguments")
		}
	}

	// Loop over every argument provided
	for ((i, arg) in args.withIndex()) {
		// If the argument is a flag
		if (arg.startsWith("-")) {
			// Checks the argument has optional arguments
			if (i == args.size - 1 || args[i + 1] in validFlags) {
				if (arg in validFlags) {
					error("Missing optional argument for $arg")
				} else {
					error("Unrecognized flag: $arg")
				}
			} 

			var optArg = args[i + 1]

			/* Check every argument different from -n have just
			 one optional argument */
			if (i + 2 != args.size && arg != "-n" && !args[i + 2].startsWith("-")) {
					error("$arg must have just one argument")
			}

			// Checks the conditions for each flag
			when (arg) {
				"-t" -> {
					if (isInt(optArg) && optArg.toInt() > 0) {
						tries = optArg.toInt()
					} else {
						error("Argument for -t must be a positive integer")
					}
				}

				"-s" -> {
                    val seqIds = arrayOf("random", "sorted", "inv",
                        "cerouno", "mitad", "repetido",)
					if (optArg in seqIds) {
						parsedArgs[0] = optArg
					} else {
						error("""Invalid optional argument for -s: $optArg.
							|Valid options are: ${seqIds.joinToString(", ")}."""
						.trimMargin("|"))
					}
				}

				"-n" -> {
					for (j in (i + 1) until args.size) {
						if (args[j] in validFlags) break

						if (!isInt(args[j]) || args[j].toInt() < 0) {
							error("Every argument for -n must be a nonnegative integer")
						} else {
							seqSizes[j] = args[j].toInt()
						}
					}
				}

				"-a" -> {
                    val algIds = arrayOf("all", "On2", "nlgn", "ln")
					if (optArg in algIds) {
						parsedArgs[1] = optArg
					} else {
						error("""Invalid optional argument for -a: $optArg.
							|Valid options are: ${algIds.joinToString(", ")}"""
						.trimMargin("|"))
					}
				}

				"-g" -> {
					if (!(optArg.endsWith(".png"))) {
						parsedArgs[2] = optArg
					} else {
						error("""Argument for -g must not end with '.png'""")
					}
				}

				else -> error("Unrecognized flag: $arg")
			}
		}
	}

	val seqSizesReturn = seqSizes.filterNotNull().toTypedArray()

	if (seqSizesReturn.size == 1 && "-g" in args) {
		error("Cannot generate a graph with just one argument for -g")
	}

	if (seqSizesReturn.size > 1 && "-g" !in args) {
		error("Must provided -g for more than two -n arguments")
	}

	seqSizesReturn.sort()

	for (i in 0 until (seqSizesReturn.size - 1)) {
		if (seqSizesReturn[i] == seqSizesReturn[i + 1]) {
			error("Every argument for -n must be distinct")
		} 
	}
	return Triple(tries, parsedArgs, seqSizesReturn)
}

/** Returns a boolean indicating is [value] in an integer */
fun isInt(value: String) =
    if (value.toIntOrNull() == null) false else true

/**
 * Checks if [A] is ordered in ascending order. Halts the 
 * program execution if not, printing an error message.
 */
fun checkSort(A: Sequence) {
    val isSorted = A.indices.all {
        it == A.size - 1 || A[it] <= A[it + 1]
    }

    if (!isSorted) {
        println("Error: sequence was not correctly  sorted")
        exitProcess(1)
    }
}

/** Returns the average of the elements of [A]. */
fun calcArrayAvg(A: Array<Long>): Double = A.sum() / A.size.toDouble()

/** Returns a pair with the smaller and the largest element of [A]. */
fun calcMinMax(A: Array<Long>): Pair<Double, Double> {
    var min = POSITIVE_INFINITY
    var max = NEGATIVE_INFINITY
    for (element in A) {
        if (element < min) {
            min = element.toDouble()
        }

        if (element > max) {
            max = element.toDouble()
        }
    }
    return Pair(min, max)
}

/**
* Generates all the necessary data to plot the
* running times of the sorting algorithms. Every
* sorting algorithm according to [algset] is
* executed [tries] times over different sequence
* sizes contained in [seqSizes] and generated 
* according to [sequence].
*
* @return A triple with the algorithm labels, the
* sizes of the input sequences, and a triple with 
* the minimun, average and maximum execution times
* of the algorithms. 
*/
fun generateData(tries: Int, sequence: String, seqSizes: Sequence, algSet: String):
    Triple<Array<String>, Sequence, Triple<Array<Double>, Array<Double>, Array<Double>>> {

    val algNames = arrayOf(
        "Bubble Sort", "Insertion Sort", "Selection Sort", "Shell Sort",
        "Merge-Insertion Sort", "Merge-Iterative Sort", "Heap Sort",
        "Randomized Quicksort", "3-way Partition Quicksort",
        "Quicksort", "Introsort", "Counting Sort", "Radix Sort",
        )

    val sortFunctions = arrayOf(
        ::bubbleSort, ::insertionSort, ::selectionSort, ::shellSort,
        ::mergeSortInsertion, ::iterativeMergeSort, ::heapSort,
        ::randomizedQuicksort, ::threeWayQuicksort,
        ::quicksort, ::introsort, ::countingSort, :: radixSort,
        )

    val start = when (algSet) {
        "nlgn" -> 4
        "ln" -> 9
        else -> 0
    }

    val end = when (algSet) {
        "On2" -> 4
        "nlgn" -> 11
        else -> algNames.size
    }

    // Declare the output arrays
    val size = seqSizes.size * (end - start)
    val labels = Array<String>(size, { "" })
    val elements = Sequence(size, { 0 })
    val minTimes = Array<Double>(size, { 0.0 })
    val maxTimes = Array<Double>(size, { 0.0 })
    val avgTimes = Array<Double>(size, { 0.0 })

    println("""
        |A summary of the running time stats will be printed
        |during the runtime to keep track of the execution.
        |The format wil be:
        |
        |Sequence size: <N> (<id>)
        |   <Algorithm>: <MinTime> <AvgTime> <MaxTime> ms
        """.trimMargin("|"))
    println()

    for ((i, N) in seqSizes.withIndex()) {

        // Generate the sequence to be sorted
        val sec = when (sequence) {
            "cerouno" -> Sequence(N, { (0..1).random() })
            "mitad" -> Sequence(N, { if (it + 1 <= ceil(N / 2.0)) it + 1 else N - it})
            "repetido" -> {
                val bound = ceil(N.toDouble().pow(2 / 3.0)).toInt()
                Sequence(N, { (0..bound).random() })
            }
            else -> Sequence(N, { (0..N).random() })
        }

        when (sequence) {
            "inv" -> sec.sortDescending()
            "sorted" -> sec.sort()
        }

        println("Sequence size: $N ($sequence)")

        // Tries the same array on every function
        for (j in 0 until (end - start)) {        
            val execTimes = Array<Long>(tries, { 0 })

            repeat (tries) { it ->
                val secCopy = sec.copyOf()

                execTimes[it] = measureTimeMillis { sortFunctions[j + start](secCopy) }
                checkSort(secCopy)
            }
            
            // Calculates the data for the output
            val index = i + j * seqSizes.size
            labels[index] = algNames[j + start]
            elements[index] = N
            val minMax = calcMinMax(execTimes)
            minTimes[index] = minMax.first
            maxTimes[index] = minMax.second
            avgTimes[index] = calcArrayAvg(execTimes)

            println("   ${labels[index]}: ${"%.2f".format(minTimes[index])} ${"%.2f".format(avgTimes[index])} ${"%.2f".format(maxTimes[index])} ms")
        }
        println()
    }
    return Triple(labels, elements, Triple(minTimes, avgTimes, maxTimes))
}

/**
 * Main program that executes a set of sorting algorithms
 * from the Sortlib.kt library on the same sequence
 * of different sizes according to the data provided by
 * the client through command line arguments.
 *
 * Opens a windows with a plot of Execution time in ms
 * vs. Number of elements and saves it im a .png file
 */
fun main(args: Array<String>) {
    // Parse the parameters provided by the command line
    val parsedArgs = parseCommands(args)

    val tries = parsedArgs.first
    val id = parsedArgs.second[0]
    val set = parsedArgs.second[1]
    val filename = parsedArgs.second[2]
    val seqSizes = parsedArgs.third

    val (labels, elements, runtimes) = generateData(tries, id, seqSizes, set)

    if (seqSizes.size != 1) {
    plotRuntime(
        "T(n) Graphics of sorting algorithms",
        ".", "$filename.png",
        "Result of the sorting algorithm in a sequence of numbers ($id)",
        "Number of elements", "Time (ms)",
        labels, elements,
        runtimes.first, runtimes.second,runtimes.third
        )
    }
}