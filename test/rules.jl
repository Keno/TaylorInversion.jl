using SymbolicUtils

@syms x::Complex y::Complex z::Complex

sqexpand = @rule (~x + ~y + ~z)/(~y) => (~x + ~z)/(~y) + 1
# sqexpand = @rule (~x + ~y + ~z)/(~y) => (~x + ~z)/(~y) + 1

# sqexpand(cos(y) + x / (x))
sqexpand((x + y + z)/y)