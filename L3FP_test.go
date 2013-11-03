package LatticeReduction

import (
	"fmt"
	"math/rand"
	"testing"
)

var SmallBasisTest = SmallBasis{
	{1, 1, 1},
	{-1, 0, 2},
	{3, 5, 6},
}

func TestSmallBasis(t *testing.T) {
	if err := SmallBasisTest.Validate(); err != nil {
		t.Error(err)
	}
}

func BenchmarkSmallBasisL3FP100x100(b *testing.B) {
	for i := 0; i < b.N; i++ {
		b.StopTimer()
		basis := make(SmallBasis, 100)
		for j := 0; j < 100; j++ {
			basis[j] = make([]int64, 100)
			for k := 0; k < 100; k++ {
				basis[j][k] = rand.Int63n(0x7FFFFFFFFFFFFFF)
			}
		}
		b.StartTimer()
		_ = basis.L3FP(0.99)
	}
}


func ExampleReduceL3FP() {
	fmt.Println(SmallBasisTest)
	fmt.Println(SmallBasisTest.L3FP(0.9))
	// Output:
	// [[1 1 1] [-1 0 2] [3 5 6]]
	// [[0 1 0] [1 0 1] [-1 0 2]]
}

func ExampleReduceL3FPDeep() {
	fmt.Println(SmallBasisTest)
	fmt.Println(SmallBasisTest.L3FPDeep(0.9))
	// Output:
	// [[1 1 1] [-1 0 2] [3 5 6]]
	// [[0 1 0] [1 0 1] [-1 0 2]]
}
