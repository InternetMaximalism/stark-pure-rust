template X() {
    signal input x;
    signal output y;
    signal x2;
    signal x3;
    var a;
    {
        a = (x*x*x-2*x*x+6)/x;
        y <-- a;
    }

    x2 <== x*x;
    x3 <== x2*x;
    x*y === x3-2*x2+6;
}

component main = X();
