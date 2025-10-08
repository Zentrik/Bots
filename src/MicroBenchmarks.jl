f(a, b) = a.*b
a = rand(Bool, 9, 12)
a1 = BitMatrix(a)
b = rand(12)

using BenchmarkTools
@benchmark f($a, $b)
@benchmark f($a1, $b)
@benchmark f($(Float64.(a)), $b)

using StaticArrays
b1 = SVector{12}(b)
a2 = SMatrix{9, 12, Bool}(a)

@benchmark f($a2, $b)

@benchmark f($a, $b1)
@benchmark f($a1, $b1)
@benchmark f($(Float64.(a)), $b1)
@benchmark f($a2, $b1)


g(a, b) = a.*b
b = rand(9)
b1 = SVector{9}(b)
a2i = SVector{9, Bool}(a[:, 2])
@benchmark g($(a[:, 2]), $b)
@benchmark g($(BitVector(a[:, 2])), $b)
@benchmark g($(Float16.(a[:, 2])), $b)
@benchmark g($(Float32.(a[:, 2])), $b)
@benchmark g($(Float64.(a[:, 2])), $b)
@benchmark g($a2i, $b)

@benchmark g($(Int16.(a[:, 2])), $b)
@benchmark g($(Int32.(a[:, 2])), $b)
@benchmark g($(Int.(a[:, 2])), $b)

@benchmark g($(a[:, 2]), $b1)
@benchmark g($(BitVector(a[:, 2])), $b1)
@benchmark g($(Float64.(a[:, 2])), $b1)
@benchmark g($a2i, $b1)

@benchmark g($(@view a[:, 2]), $b)

test = Float32.(a)
@benchmark g($(@view test[:, 2]), $b)

function test2()
    test = Float32.(a)
    @benchmark g($@view(test[:, 2]), $b)
end