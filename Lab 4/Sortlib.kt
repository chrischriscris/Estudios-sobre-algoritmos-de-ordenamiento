import kotlin.math.*
typealias Sequence = Array<Int>

/** Returns the greatest integer less or equal than lg(n). */
fun floorLg(n: Int): Int = floor(log2(n.toDouble())).toInt()

/** Returns the median of [a], [b] and [c]. */
fun medianOf3(a: Int, b: Int, c: Int): Int = arrayOf(a, b, c).sorted()[1]

/** 
 * Prints the elements of [A] in the range [start, end).
 * For debugging purposes.
 */
fun printSequence(A: Sequence, start: Int = 0, end: Int = A.size) {
    for (i in start until end) {
        print("${A[i]} ")
    }
    println()
}

/** Exchanges [A[i]] and [A[j]] */
fun swap(A: Sequence, i: Int, j: Int) {
    val temp = A[i]
    A[i] = A[j]
    A[j] = temp
}

/** Sorts [A] in increasing order using Bubble Sort. */
fun bubbleSort(A: Sequence) {
    val n = A.size

    for (i in 0 until (n - 1)) {
        for (j in (n - 1) downTo (i + 1)) {
            if (A[j - 1] > A[j]) {
                swap(A, j - 1, j)
            }
        }
    }
}

/** Sorts [A] in range [l, r] in increasing order using Insertion Sort. */
fun insertionSortRange(A: Sequence, l: Int, r: Int) {
    for (i in (l + 1)..r) {
        var j = i
        while (j != l && A[j] < A[j - 1]) {
            swap(A, j, j - 1)
            j--
        }
    }
}

fun insertionSort(A: Sequence) = insertionSortRange(A, 0, A.size - 1)

/** Sorts [A] in increasing order using Selection Sort. */
fun selectionSort(A: Sequence) {
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
        swap(A, i, lowindex)
    }
}

/** Sorts [A] in increasing order using Shell Sort. */
fun shellSort(A: Sequence) {
    val n = A.size

    var incr = n / 2
    while (incr > 0) {
        for (i in incr until n) {
            var j = i - incr

            while (j >= 0) {
                if (A[j] > A[j + incr]) {
                    swap(A, j, j + incr)
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
 * m and n, respectively, into one ordered array [A], of size m+n.
 */
fun merge(U: Sequence, V: Sequence, A: Sequence) {
    var (i, j) = Pair(0, 0)

    val (m, n) = Pair(U.size, V.size)

    for (k in 0 until (m + n)) {
        if (i != m && (j == n || U[i] < V[j])) {
            A[k] = U[i]
            i++
        } else {
            A[k] = V[j]
            j++
        }
    }
}

/** Sorts [A] in increasing order using Merge-Insertion Sort. */
fun mergeSortInsertion(A: Sequence) {
    val n = A.size

    if (n <= 10) {
        insertionSort(A)
    } else {
        val U = A.sliceArray(0 until (n / 2))
        val V = A.sliceArray((n / 2) until n)
        mergeSortInsertion(U)
        mergeSortInsertion(V)
        merge(U, V, A)
    }
}

/** Sorts [A] in increasing order using an iterative Sort. */
fun iterativeMergeSort(A: Sequence) {
    val n = A.size
    var k = 1

    while (k < n) {
        var (a, b, c) = Triple(0, k, minOf(2 * k, n))

        /* Orders and merges subsequent subarrays
        of T with size of (at most) 2^(k-1) */
        while (b < n) {
            var (p, q, r) = Triple(a, b, 0)
            var Z = Sequence(c - a, { 0 })

            // Merges A[a..b) with A[b..c) in Z[0..c-a)
            while (p != b && q != c) {
                if (A[p] <= A[q]) {
                    Z[r] = A[p]
                    r++
                    p++
                } else {
                    Z[r] = A[q]
                    r++
                    q++
                }
            }

            // Terminates to empty A into Z
            while (p != b) {
                Z[r] = A[p]
                r++
                p++
            }
            while (q != c) {
                Z[r] = A[q]
                r++
                q++
            } 

            // Copies Z[0..a-c) into A[a..c)
            r = a 
            while (r != c) {
                A[r] = Z[r - a]
                r++
            }

            a += 2 * k
            b += 2 * k
            c = minOf(c + 2 * k, n)
        }
        k *= 2
    }
}

/**
 * Maintains a Max-Heap properties of a Heap in [A] in [start, end)
 * out of parent node in index [i].
 */
fun maxHeapify(A: Sequence, i: Int, start: Int, end: Int) {
    var (l, r) = Pair(2 * i + 1 - start, 2 * i + 2 - start)
    var largest: Int

    largest = if (l < end && A[l] > A[i]) l else i
    largest = if (r < end && A[r] > A[largest]) r else largest 

    if (largest != i) {
        swap(A, i, largest)
        maxHeapify(A, largest, start, end)
    }
}

/** Builds a Max-Heap out of all the elements in [A] in range [l,r). */
fun buildMaxHeap(A: Sequence, l: Int, r: Int) {
    val n = r - l
    for (i in (n / 2 - 1 + l)  downTo l) {
        maxHeapify(A, i, l, r)
    }
}

/** Sorts [A] in [l, r) in increasing order using Heap Sort. */
fun heapSortRange(A: Sequence, l: Int, r: Int) {
    buildMaxHeap(A, l, r)
    var end = r
    for (i in (r - 1) downTo (l + 1)) {
        swap(A, i, l)
        end--
        maxHeapify(A, l, l, end)
    }
}

fun heapSort(A: Sequence) = heapSortRange(A, 0, A.size)

/**
 * Rearranges the elements in [A] in range [p, r] and returns
 * an l such that A[k] <= A[r] for all p <= k < l and
 * A[k] > A[r] for all l < k <= r.
 */
fun partition(A: Sequence, p: Int, r: Int): Int {
    val x = A[r]
    var i = p - 1

    for (j in p until r) {
    	if (A[j] <= x) {
            i++
            swap(A, i, j)
        }
    }
    swap(A, i + 1, r)
    return i + 1
}

/** Sorts [A] in range [p, r] in increasing order using Quicksort */
fun quicksortRange(A: Sequence, p: Int, r: Int) {
    if (p < r) {
    	val q = partition(A, p, r)
    	quicksortRange(A, p, q - 1)
    	quicksortRange(A, q + 1, r)
    }
}

/** Sorts [A] in increasing order using Quicksort */
fun quicksort(A: Sequence) = quicksortRange(A, 0, A.size - 1) 

/**
 * Rearranges the elements in [A] in range [p, r] and returns
 * an l such that A[k] <= A[c] for all p <= k < l and
 * A[k] > A[c] for all l < k <= r, where c is a random
 * number between p and r, both inclusive.
 */
fun randomizedPartition(A: Sequence, p: Int, r: Int): Int {
    val i = (p..r).random()
    swap(A, r, i)
    return partition(A, p, r)
}

/**
 * Sorts [A] in range [p, r] in increasing order using a variant
 * of Quicksort with random partitioning.
 */
fun randomizedQuicksortRange(A: Sequence, p: Int, r: Int) {
    if (p < r) {
        val q = randomizedPartition(A, p, r)
        randomizedQuicksortRange(A, p, q - 1)
        randomizedQuicksortRange(A, q + 1, r)
    }
}

/** 
 * Sorts [A] in increasing order using a variant of Quicksort with
 * random partitioning.
 */
fun randomizedQuicksort(A: Sequence) = randomizedQuicksortRange(A, 0, A.size - 1)

/**
 * Sorts [A] in range [p, r] in increasing order using a variant of
 * Quicksort with 3-way partitioning.
 */
fun threeWayQuicksortRange(A: Sequence, l: Int, r: Int) {

    if (r <= l) return

    var (i, j) = Pair(l - 1, r)
    var (p, q) = Pair(l - 1, r)
    val v = A[r]

    while (true) {

        /* Find the index i of the first element from
        A[l] to A[r] to be greater than the pivot and
        the first one from r to l to be smaller */
        while (A[++i] < v) {}
        while (A[--j] > v) if (j == l) break
        
        if (i >= j) break
        swap(A, i, j)

        /*If any of both is equal to the pivot,
        then swap it with a corner */
        if (A[i] == v) {
            p++
            swap(A, p, i)
        }
        if (A[j] == v) {
            q--
            swap(A, j, q)
        }

    }

    // Make A[i] the pivot
    swap(A, i, r)
    j = i - 1
    i++

    // Move the corners to the center
    for (k in l until p){
        swap(A, k, j)
        j--
    }
    for (k in (r - 1) downTo (q + 1)) {
        swap(A, i, k)
        i++
    }

    threeWayQuicksortRange(A, l, j)
    threeWayQuicksortRange(A, i, r)
}

/**
 * Sorts [A] in increasing order using a variant of Quicksort
 * with 3-way partitioning.
 */
fun threeWayQuicksort(A: Sequence) = threeWayQuicksortRange(A, 0, A.size - 1)

/**
 * Rearranges the elements in [A] in range [p..r) and returns
 * an l such that A[k] <= x for all p <= k <= l and
 * A[k] >= x for all l < k < r.
 */
fun hoarePartition(A: Sequence, p: Int, r: Int, x: Int): Int {
    var (i, j) = Pair(p - 1, r)

    while (true) {
        do { j-- } while (A[j] > x)
        do { i++ } while (A[i] < x)

        if (i < j) swap(A, i, j) else return i
    }
}

/**
 * Sorts [A] in range [f, b) in increasing order using Introspective Sort.
 * Uses Quicksort pivoting on the median of A[f], A[r] and A[(f+r) / 2] and
 * limits the depth of recursive calls to [depthLimit]. Orders the
 * remaining portions when [depthLimit] is exceeded using Heap Sort.
 */
fun introsortLoop(A: Sequence, f: Int, b: Int, depthLimit: Int) {
    var r = b
    val sizeThreshold = 32

    while (r - f > sizeThreshold) {
        if (depthLimit == 0) {
            heapSortRange(A, f, r)
            return
        }

        val x = medianOf3(A[f], A[f + ((r - f) / 2)], A[r - 1])
        val p = hoarePartition(A, f, r, x)
        introsortLoop(A, p, r, depthLimit - 1)
        r = p
    }
}

/** Sorts [A] in range [f, b) increasing order using Introspective Sort */
fun introsortRange(A: Sequence, f: Int, b: Int) {
    introsortLoop(A, f, b, 2 * floorLg(b - f))
    insertionSortRange(A, f, b - 1)
}

/** Sorts [A] in increasing order using Introspective Sort */
fun introsort(A: Sequence) = introsortRange(A, 0, A.size)