package LatticeReduction

import (
	"fmt"
	"math/big"
)

// Lattice basis for big.Int sized values
type BigBasis [][]*big.Int

var (
	bigOne = big.NewInt(1)
)

func bigIntToFloat(in *big.Int) float64 {
	rat := new(big.Rat)
	rat.SetFrac(in, bigOne)
	f, _ := rat.Float64()
	return f
}

func (b BigBasis) Copy()Basis {
	o:=make(BigBasis, len(b))
	for i, x := range b {
		o[i] = make([]*big.Int,len(x))
		for j,y:=range x{
			o[i][j] = new(big.Int)
			o[i][j].Set(y)
		}
	}
	return o
}

func (b BigBasis) ColumnSwap(column1,column2 int){
	b[column1], b[column2] = b[column2], b[column1]
}

func (b BigBasis) Rank() int{
	return len(b)
}

func (b BigBasis) Dimension() int{
	return len(b[0])
}

func (b BigBasis) FDot(column1, column2 int)(r float64){
	var intermediate = new(big.Int)
	for i, x := range b[column1] {
		intermediate.Mul(x, b[column2][i])
		r += bigIntToFloat(intermediate)
	}
	return
}

func (b BigBasis) FColumnReduce(column1, column2 int, mu float64){
	var (
		intermediate = new(big.Int)
		_mu          = big.NewInt(int64(mu))
	)

	for i := range b[column1] {
		intermediate.Mul(_mu, b[column2][i])
		b[column1][i].Sub(b[column1][i], intermediate)
	}
}

func (b BigBasis) FGet(column,row int)float64{
	return bigIntToFloat(b[column][row])
}

func (b BigBasis) String() (s string) {
	for _, r := range b {
		s += fmt.Sprintln(r)
	}
	return s
}
