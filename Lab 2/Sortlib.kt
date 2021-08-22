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

/**
 * Merges two ordered subarrays [U] and [V] of a given size
 * m and n, respectively, into one ordered array [T], of size m+n.
 */
fun merge(U: Array<Int>, V: Array<Int>, T: Array<Int>) {
    var i = 0
    var j = 0

    val m = U.size
    val n = V.size

    for (k in 0 until (m + n)) {
        if (i != m && (j == n || U[i] < V[j])) {
            T[k] = U[i]
            i++
        } else {
            T[k] = V[j]
            j++
        }
    }
}

/** Sorts [A] in increasing order using Merge-Insertion Sort. */
fun mergeSortInsertion(T: Array<Int>) {
    val n = T.size

    if (n <= 10) {
        insertionSort(T)
    } else {
        val U = T.sliceArray(0 until (n / 2))
        val V = T.sliceArray((n / 2) until n)
        mergeSortInsertion(U)
        mergeSortInsertion(V)
        merge(U, V, T)
    }
}

/** Sorts [A] in increasing order using Merge-Iterative Sort. */
fun mergeSortIterativo(T: Array<Int>) {
    val n = T.size

    var k = 1
    while (k < n) {
        var a = 0
        var b = k
        var c = minOf(2 * k, n)

        /* Orders and merges subsequent subarrays
        of T with size of (at most) 2^(k-1) */
        while (b < n) {
            var p = a
            var q = b
            var r = 0
            var Z = Array<Int>(c - a, { 0 })

            // Merges T[a..b) with T[b..c) in Z[0..c-a)
            while (p != b && q != c) {
                if (T[p] <= T[q]) {
                    Z[r] = T[p]
                    r++
                    p++
                } else {
                    Z[r] = T[q]
                    r++
                    q++
                }
            }

            // Terminates to empty T into Z
            while (p != b) {
                Z[r] = T[p]
                r++
                p++
            }
            while (q != c) {
                Z[r] = T[q]
                r++
                q++
            } 

            // Copies Z[0..a-c) into T[a..c)
            r = a 
            while (r != c) {
                T[r] = Z[r - a]
                r++
            }

            a += 2 * k
            b += 2 * k
            c = minOf(c + 2 * k, n)
        }
        k *= 2
    }
}