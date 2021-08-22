import kotlin.system.measureTimeMillis
import kotlin.math.pow
import kotlin.math.sqrt
import kotlin.system.exitProcess

typealias Sequence = Array<Int>

/**
 * Returns the execution time in milliseconds
 * of the [sortFunction] applied to the given [A].
 */
fun measureExecutionTime(
    sortFunction: (Sequence) -> Unit,
    A: Sequence,
    ): Long = measureTimeMillis {
        sortFunction(A)
    }

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

/**
 * Returns a pair with the average and the standar deviation
 * of the elements of [A].
 */
fun calcArrayStats(A: Array<Long>): Pair<Double, Double> {
    val n = A.size
    val avg = A.sum() / n.toDouble()
    var stdDev = 0.0 
    
    for (x in A) {
        var distance = x - avg
        stdDev += distance.pow(2)
    }

    stdDev /= n
    stdDev = sqrt(stdDev)

    return Pair(avg, stdDev)
}

/**
 * Tries Insertion, Selection, Shell and Bubble Sort over
 * a randomly generated sequence a certain amount of times
 * indicated through [args].
 *
 * Keeps track of execution times and prints it to the standard
 * output. If there is more than one try, the output will
 * reflect the average execution time and the standard
 * deviation.
 */
fun main(args: Array<String>) {
    // Parse the parameters provided by the command line
    val tries = args[0].toInt()
    val id = args[1]
    val N = args[2].toInt()

    // Generate the sequence to be sorted
    val sec = Sequence(N, { (0..N).random() })
    when (id) {
        "inv" -> sec.sortDescending()
        "sorted" -> sec.sort()
    }

    val sortFunctions = arrayOf(
        ::insertionSort,
        ::selectionSort,
        ::shellSort,
        ::bubbleSort,
        )

    val names = arrayOf("Insertion", "Selection", "Shell", "Bubble",)

    println("""Execution time metrics:
        | (For $tries tries on the same sequence)
        | (Size of the sequence: $N)
        | (Type of sequence: $id)
        """.trimMargin("|"))
    println()

    /* Loop over all the functions in the array
    and tries the sort algorithm */
    for ((sortFunction, name) in sortFunctions zip names) {
        if (tries == 1) {
            var secCopy = sec.copyOf()
            var time = measureExecutionTime(sortFunction, secCopy)
            checkSort(secCopy)
            println("$name Sort: $time ms")
        } else {
            var execTimes = Array<Long>(tries, { 0 })

            println("$name sort:\n")
            repeat (tries) { it ->
                var secCopy = sec.copyOf()
                execTimes[it] = measureExecutionTime(sortFunction, secCopy)
                println("  ${it + 1}º try: ${execTimes[it]} ms")
                checkSort(secCopy)
            }

            if (tries != 0) {
                val stats = calcArrayStats(execTimes)

                println()
                println("""
                    |  Average running time: ${"%.2f".format(stats.first)} ms
                    |  Standard deviation: ${"%.2f".format(stats.second)} ms
                    """.trimMargin("|"))
                println()
            }
        }
    }
}