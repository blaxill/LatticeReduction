package LatticeReduction

type Int64Vector[]int64
// type SmallBasis []Int64Vector

func (lhs Int64Vector) Len() int{
	return len(lhs)
}

// func (lhs Int64Vector) FDot(rhs Float64Vector) float64 {
// 	for i, x := range lhs {
// 		r += float64(x * rhs[i])
// 	}
// 	return
// }

func (lhs Int64Vector) Dot(rhs Vector) (r float64) {
	_rhs, ok := rhs.(Int64Vector)
	if !ok{
		panic("Something went wrong.")
	}
	for i, x := range lhs {
		r += float64(x * _rhs[i])
	}
	return
}

func (lhs Int64Vector) FReduceInplace(mu float64, rhs Vector) {
	_rhs, ok := rhs.(Int64Vector)
	if !ok{
		panic("Something went wrong.")
	}
	for i := range lhs {
		lhs[i] -= int64(mu * float64(_rhs[i]))
	}
}

func (lhs Int64Vector) Copy() Vector {
	r:=make(Int64Vector,len(lhs))
	for i, x := range lhs {
		r[i] = x
	}
	return r
}

func (lhs Int64Vector) FAt(i int) float64 {
	return float64(lhs[i])
}
