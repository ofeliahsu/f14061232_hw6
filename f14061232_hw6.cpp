// f14061232_hw6.cpp : 此檔案包含 'main' 函式。程式會於該處開始執行及結束執行。
//

#include <iostream>
#include <cstdlib>
#include <vector>
#include<list>
#include <algorithm>
#include <stdio.h>      /* printf */
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <math.h>       /* sqrt */
#include <cstdio>
#include <ctime>

using namespace::std;

void reverse(int arr[], int n)
{
    reverse(arr, arr + n);
}


int GenRanInt(int a_min, int a_max) {
    return (a_min + rand() % (a_max - a_min + 1));
}

int is_A_desc(int* A, int n) {
    for (int i = 0; i < n; i++)
    {
        int j=1;
        if (A[i] >= A[j] && i < j)
        {
            cout << "order OK!";
            return 1;
        }
        else {
            cout << "!!WRONGLY SORTED!!";
            return 0;
        }
    }
}

void printout_arrary(int* A, int n) {
    for (int i = 0; i < n; ++i)
    {
            cout << A[i] << ' ';
    }
    cout << endl;
}

void duplicate_array(int* A, int* B, int n) {

        for (int i = 0; i < n; i++)
        {
            B[i] = A[i];
        }

}



void ins_sort(int* A, int n)
{
    int temp = 0;
    int next = 0;

    for (int i = 0; i < n; i++)
    {
        temp = A[i];//取出一個key
        next = i;//已排序中最後一個數值的位置

        while ((next > 0) && (A[next - 1] > temp))//已排序中最後一個數值的位置>0且其數值比key大時
        {
            A[next] = A[next - 1];//將key前一個數值複製給key的位置
            next = next - 1;//key繼續往前比較大小
        }

        A[next] = temp;//將key插入排序中其值<=key的元素後面

    }
}

void merge(int* B, int l, int m, int r)
{
    int n1 = m - l + 1;
    int n2 = r - m;

    // Create temp arrays
    int* L = new int[n1];
    int* R = new int[n2];

    // Copy data to temp arrays L[] and R[]
    for (int i = 0; i < n1; i++)
        L[i] = B[l + i];

    for (int j = 0; j < n2; j++)
        R[j] = B[m + 1 + j];

    // Merge the temp arrays back into arr[l..r]

    // Initial index of first subarray
    int i = 0;

    // Initial index of second subarray
    int j = 0;

    // Initial index of merged subarray
    int k = l;

    while (i < n1 && j < n2)
    {
        if (L[i] <= R[j])
        {
            B[k] = L[i];
            i++;
        }
        else
        {
            B[k] = R[j];
            j++;
        }
        k++;
    }

    while (i < n1)
    {
        B[k] = L[i];
        i++;
        k++;
    }

    while (j < n2)
    {
        B[k] = R[j];
        j++;
        k++;
    }

}


void mg_sort(int* B, int l, int r)
{
    if (l < r)
    {

        int m = l + (r - l) / 2;
        mg_sort(B, l, m);
        mg_sort(B, m + 1, r);
        merge(B, l, m, r);
    }

}

void heapify(int* A, int n, int i)
{
    int largest = i; // Initialize largest as root
    int l = 2 * i + 1; // left = 2*i + 1
    int r = 2 * i + 2; // right = 2*i + 2

    // If left child is larger than root
    if (l < n && A[l] > A[largest])
        largest = l;

    // If right child is larger than largest so far
    if (r < n && A[r] > A[largest])
        largest = r;

    // If largest is not root
    if (largest != i) {
        swap(A[i], A[largest]);

        // Recursively heapify the affected sub-tree
        heapify(A, n, largest);
    }

}



void heap_sort(int* A, int n) {
    // Build heap (rearrange array)
    for (int i = n / 2 - 1; i >= 0; i--)
        heapify(A, n, i);

    // One by one extract an element from heap
    for (int i = n - 1; i >= 0; i--) {
        // Move current root to end
        swap(A[0], A[i]);

        // call max heapify on the reduced heap
        heapify(A, i, 0);
    }


}

int Partition(int* B, int front, int end) {
    int pivot = B[end];
    int i = front - 1;
    for (int j = front; j < end; j++) {
        if (B[j] < pivot) {
            i++;
            swap(B[i], B[j]);
        }
    }
    i++;
    swap(B[i], B[end]);
    return i;
}

void quick_sort(int* B, int front, int end) {

    if (front < end) {
        int pivot = Partition(B, front, end);
        quick_sort(B, front, pivot - 1);
        quick_sort(B, pivot + 1, end);
    }

}

void countingSort(int array[], int size) {

    int output[10];
    int count[10];
    int max = array[0];

    // Find the largest element of the array
    for (int i = 1; i < size; i++) {
        if (array[i] > max)
            max = array[i];
    }

    // Initialize count array with all zeros.
    for (int i = 0; i <= max; ++i) {
        count[i] = 0;
    }

    // Store the count of each element
    for (int i = 0; i < size; i++) {
        count[array[i]]++;
    }

    // Store the cummulative count of each array
    for (int i = 1; i <= max; i++) {
        count[i] += count[i - 1];
    }

    // Find the index of each element of the original array in count array, and
    // place the elements in output array
    for (int i = size - 1; i >= 0; i--) {
        output[count[array[i]] - 1] = array[i];
        count[array[i]]--;
    }

    // Copy the sorted elements into original array
    for (int i = 0; i < size; i++) {
        array[i] = output[i];
    }

}

void bucketSort(int array[], int n)
{
    vector<int>b[n];

    for (int i = 0; i < n; i++)
    {
        int x = n * array[i];
        b[x].push_back(array[i]);
    }

    for (int i = 0; i < n; i++)
        sort(b[i].begin(), b[i].end());

    int index = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < b[i].size(); j++)
            array[index++] = b[i][j];

}

void radixSort(int* array, int n) {
    int i, j, m, p = 1, index, temp, count = 0;
    int max = array[0];

    for (int i = 1; i < n; i++) {
        if (array[i] > max)
            max = array[i];
    }
    list<int> pocket[10];      //radix of decimal number is 10
    for (i = 0; i < max; i++) {
        m = pow(10, i + 1);
        p = pow(10, i);
        for (j = 0; j < n; j++) {
            temp = array[j] % m;
            index = temp / p;      //find index for pocket array
            pocket[index].push_back(array[j]);
        }
        count = 0;
        for (j = 0; j < 10; j++) {
            //delete from linked lists and store to array
            while (!pocket[j].empty()) {
                array[count] = *(pocket[j].begin());
                pocket[j].erase(pocket[j].begin());
                count++;
            }
        }
    }

}



int main()
{
    int i, j, k;
    int is_p, n, a_min, a_max, rs;
    float t_is, t_ms, t_hs, t_qs, t_cs, t_rs, t_bs;
    std::clock_t start;
    double duration;


    cout << "Input [is_p,n,a_min,a_max,rs]= ";
    cin >> is_p >> n >> a_min >> a_max >> rs;
    srand(rs);
    int* A= new int[n];
    int* B = new int[n];

    for (i = 0; i < n; ++i) {
        A[i] = GenRanInt(a_min, a_max);
    }

    duplicate_array(A, B, n);
start = std::clock();
    ins_sort(A, n);
t_is = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    reverse(A, n);

    is_A_desc(A, n);
    while (is_A_desc(A, n) == 0) {
        cout << "!!INSERTION SORT WRONG!!";
        break;
    }
    if (is_p == 1) {
        printout_arrary(A, n);
    }


    duplicate_array(B, A, n);
start = std::clock();
    mg_sort(A, 0, n-1);
 t_ms= ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

      reverse(A, n);

    is_A_desc(A, n);
    while (is_A_desc(A, n) == 0) {
        cout << "!!MERGE SORT WRONG!!";
        break;
    }
    if (is_p == 1) {
        printout_arrary(A, n);
    }


    duplicate_array(A, B, n);
start = std::clock();
    heap_sort(A, n);
 t_hs = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;


        reverse(A, n);

    is_A_desc(A, n);
    while (is_A_desc(A, n) == 0) {
        cout << "!!HEAP SORT WRONG!!";
        break;
    }
    if (is_p == 1) {
        printout_arrary(A, n);
    }


    duplicate_array(B, A, n);
start = std::clock();
    quick_sort(A, 0, n - 1);
  t_qs = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        reverse(A, n);

    is_A_desc(A, n);
    while (is_A_desc(A, n) == 0) {
        cout << "!!QUICK SORT WRONG!!";
        break;
    }
    if (is_p == 1) {
        printout_arrary(A, n);
    }


    duplicate_array(B, A, n);
start = std::clock();
    countingSort(A, n);
  t_cs = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

        reverse(A, n);

        is_A_desc(A, n);
    while (is_A_desc(A, n) == 0) {
        cout << "!!COUNTING SORT WRONG!!";
        break;
    }
    if (is_p == 1) {
        printout_arrary(A, n);
    }

start = std::clock();
    radixSort(A, n);
 t_rs = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        reverse(A, n);

    is_A_desc(A, n);
    while (is_A_desc(A, n) == 0) {
        cout << "!!RADIX SORT WRONG!!";
        break;
    }
    if (is_p == 1) {
        printout_arrary(A, n);
    }

start = std::clock();
    bucketSort(A, n);
 t_bs = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    reverse(A, n);

    is_A_desc(A, n);
    while (is_A_desc(A, n) == 0) {
        cout << "!!BUCKET SORT WRONG!!";
        break;
    }
    if (is_p == 1) {
        printout_arrary(A, n);
    }

cout << n<<' '<<a_min << ' ' << a_max << ' ' << rs << ' ' << t_is << ' ' << t_ms << ' ' << t_hs << ' ' << t_qs << ' ' << t_cs << ' ' << t_rs << ' ' << t_bs;

    return 0;
}
