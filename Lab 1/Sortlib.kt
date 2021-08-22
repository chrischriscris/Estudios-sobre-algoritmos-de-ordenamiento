/** Sorts [A] in increasing order using Bubble Sort. */
fun bubbleSort(A: Array<Int>) {
    val n = A.size

    for (i in 0 until (n - 1)) {
        for (j in (n - 1) downTo (i + 1)) {
            if (A[j - 1] > A[j]) {
                // Swap A[j-1] and A[j]
                var temp = A[j - 1]
                A[j - 1] = A[j]
                A[j] = temp
            }
        }
    }
}

/** Sorts [A] in increasing order using Insertion Sort. */
fun insertionSort(A: Array<Int>) {
    val n = A.size

    for (i in 1 until n) {
        var j = i
        while (j != 0 && A[j] < A[j - 1]) {
            // Swap A[j] and [j-1]
            var temp = A[j]
            A[j] = A[j - 1]
            A[j - 1] = temp
            j--
        }
    }
}

/** Sorts [A] in increasing order using Selection Sort. */
fun selectionSort(A: Array<Int>) {
    val n = A.size
    var lowkey: Int
    var lowindex: Int

    for (i in 0 until (n - 1)) {
        lowindex = i
        lowkey = A[i]

        for (j in (i + 1) until n) {
            // Compares each element with lowkey
            if (A[j] < lowkey) {
                lowkey = A[j]
                lowindex = j
            }
        }
        // Swap A[i] and A[lowindex]
        var temp = A[i]
        A[i] = A[lowindex]
        A[lowindex] = temp
    }
}

/** Sorts [A] in increasing order using Shell Sort. */
fun shellSort(A: Array<Int>) {
    val n = A.size

    var incr = n / 2
    while (incr > 0) {
        for (i in incr until n) {
            var j = i - incr

            while (j >= 0) {
                if (A[j] > A[j + incr]) {
                    // Swaps A[j] and A[j + incr]
                    var temp = A[j]
                    A[j] = A[j + incr]
                    A[j + incr] = temp
                    j -= incr
                } else {
                    break
                }
            }
        }
        incr /= 2
    }
}