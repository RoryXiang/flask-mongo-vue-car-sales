var math = require('mathjs');

DifferentialEquations = {
    /*
     * Approach to solving the initial value problem using the
     * Runge-Kutta formula of order 4.
     *
     * @param {function} f - function that we want to find the solution
     * @param {Number|BigNumber} a - right limit
     * @param {Number|BigNumber} b - left limit
     * @param {Number|BigNumber} ya - initial condition y(a)
     * @param {Number} m - number of steps
     * @param {Object} result - return 2 elements, T a abscissa vector and
     *                          Y a ordinate vector.
     * */
    rungeKuttaOrder4 : function(f, a, b, ya, m) {
        var h = math.divide(math.subtract(b, a), m);
        var Y = math.zeros(m + 1);
        //var T = math.zeros(m + 1);
        var T = math.range(a, b, h);
        Y = math.subset(Y, math.index(0), ya);

        for(var j = 0; j < m; ++j) {
            var Tj = math.subset(T, math.index(j));
            var Yj = math.subset(Y, math.index(j));

            var k1 = math.multiply(h, f(Tj, Yj));
            var k2 = math.multiply(h, f(math.add(Tj, math.divide(h, 2)), math.add(Yj, math.divide(k1, 2))));
            var k3 = math.multiply(h, f(math.add(Tj, math.divide(h, 2)), math.add(Yj, math.divide(k2, 2))));
            var k4 = math.multiply(h, f(math.add(Tj, h), math.add(Yj, k3)));

            var sum = math.add(k1, math.add(math.multiply(2, k2), math.add(math.multiply(2, k3), k4)));
            var r = math.add(Yj, math.divide(sum, 6));
            Y = math.subset(Y, math.index(j + 1), r);
        }
        var result = {T: T, Y: Y};
        return result;
    },

    /*
     * Approach to solving the initial value problem using the
     * Runge-Kutta formula of order 2.
     *
     * @param {function} f - function that we want to find the solution
     * @param {Number|BigNumber} a - right limit
     * @param {Number|BigNumber} b - left limit
     * @param {Number|BigNumber} ya - initial condition y(a)
     * @param {Number} m - number of steps
     * @param {Object} result - return 2 elements, T a abscissa vector and
     *                          Y a ordinate vector.
     * */

    rungeKuttaOrder2 : function(f, a, b, ya, m) {
        var h = math.divide(math.subtract(b, a), m);
        var Y = math.zeros(m + 1);
        //var T = math.zeros(m + 1);
        var T = math.range(a, b, h);
        Y = math.subset(Y, math.index(0), ya);

        for(var j = 0; j < m; ++j) {
            var Tj = math.subset(T, math.index(j));
            var Tj1 = math.subset(T, math.index(j + 1));
            var Yj = math.subset(Y, math.index(j));

            var k1 = f(Tj, Yj);
            var k2 = f(Tj1, math.add(Yj, math.multiply(h, k1)));

            var r = math.add(Yj, math.multiply(math.divide(h, 2), math.add(k1, k2)));
            Y = math.subset(Y, math.index(j + 1), r);
        }

        var result = {T: T, Y: Y};
        return result;
    }
};
