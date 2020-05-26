var math = require("mathjs");
require("util");

Optimization = {

    /*
    * Function that transform a "Matrix" vector into a Vector
    *
    * @param {Matrix} vect - Vector that mathjs suppose as a matrix.
    * @param {Vector} normalVector - Vector transform.
    *
    * */
    _toNormalVector : function(vect) {
        var normalVector = [];
        var size = math.subset(math.size(vect), math.index(1));
        for (var i = 0; i < size; ++i) {
            normalVector.push(math.subset(vect, math.index(0, i)));
        }
        return normalVector;
    },

    /*
    * Function that obtains the local minimum vector of a function with a vertices of a Simplex.
    *
    * * Note: function based on the book Numerical Methods by John H. Mathews and Kurtis D. Fink,
     *       Spanish, Chap 8, pag 450.
    *
    * @param {Function} f - Function to get the local minimum, this function require one argument and should be a vector.
    * @param {Vector} vertices - Vertices of the Simplex shape.
    * @param {Number} minimumIterations - Minimum number of iterations we want.
    * @param {Number} maximumIterations - Maximum number of iterations we want.
    * @param {Number|BigNumber} epsilon - Local minimum tolerance.
    * @return {Object} result - return an object that contains a minimum vector (.vector) of the function and the evaluation of the function in this vector (.eval).
    *
    * */
    This : function (f, vertices, minimumIterations, maximumIterations, epsilon) {
        var size = math.subset(math.size(vertices), math.index(1));
        var y = [];
        for(var i = 0; i < size + 1; ++i) {
            var vector = vertices[i];
            y.push(f(vector));
        }

        var low = Util.getIndexRow(y, math.min(y));
        var high = Util.getIndexRow(y, math.max(y));
        var li = high;
        var ho = low;

        for(var i = 0; i <= size; ++i) {
            if(!math.equal(i, low) && !math.equal(i, high) && math.smallerEq(y[i], y[li])) {
                li = i;
            }
            if(!math.equal(i, low) && !math.equal(i, high) && math.largerEq(y[i], y[ho])) {
                ho = i;
            }
        }

        var counter = 0;

        while( (math.larger(y[high], math.add(y[low], epsilon)) && counter < maximumIterations) || counter < minimumIterations ) {

            var S = math.zeros(1, size);
            for(var i = 0; i < size + 1; ++i) {
                S = math.add(S, math.subset(vertices, math.index(i, [0, size])));
            }

            var hiVector = math.subset(vertices, math.index(high, [0, size]));
            var M = math.divide(math.subtract(S, hiVector), size);
            var R = math.subtract(math.multiply(2, M), hiVector);
            var yR = f(this._toNormalVector(R));

            if(math.smaller(yR, y[ho])) {
                if(math.smaller(y[li], yR)) {
                    vertices = math.subset(vertices, math.index(high, [0, size]), R);
                    y[high] = yR;
                } else {
                    var E = math.subtract(math.multiply(2, R), M);
                    yE = f(this._toNormalVector(E));

                    if(math.smaller(yE, y[li])) {
                        vertices = math.subset(vertices, math.index(high, [0, size]), E);
                        y[high] = yE;
                    } else {
                        vertices = math.subset(vertices, math.index(high, [0, size]), R);
                        y[high] = yR;
                    }
                }
            } else {
                if(math.smaller(yR, y[high])) {
                    vertices = math.subset(vertices, math.index(high, [0, size]), R);
                    y[high] = yR;
                }

                var actual = math.subset(vertices, math.index(high, [0, size]));
                var C = math.divide(math.add(actual, M), 2);
                var yC = f(this._toNormalVector(C));

                var C2 = math.divide(math.add(M, R), 2);
                var yC2 = f(this._toNormalVector(C2));

                if(math.smaller(yC2, yC)) {
                    C = C2;
                    yC = yC2;
                }
                if(math.smaller(yC, y[high])) {
                    vertices = math.subset(vertices, math.index(high, [0, size]), C);
                    y[high] = yC;
                } else {
                    for(var i = 0; i <= size; ++i) {
                        if(i != low) {
                            var vectorj = math.subset(vertices, math.index(j, [0, size]));
                            var vectorlo = math.subset(vertices, math.index(low, [0, size]));

                            vertices = math.subset(vertices, math.index(j, [0, size]), math.divide(math.add(vectorj, vectorlo),2));

                            var Z = math.subset(vertices, math.index(j, [0, size]));
                            y[i] = f(this._toNormalVector(Z));
                        }
                    }
                }
            }

            low = Util.getIndexRow(y, math.min(y));
            high = Util.getIndexRow(y, math.max(y));
            li = high;
            ho = low;

            for(var i = 0; i <= size; ++i) {
                if(!math.equal(i, low) && !math.equal(i, high) && math.smallerEq(y[i], y[li])) {
                    li = i;
                }
                if(!math.equal(i, low) && !math.equal(i, high) && math.largerEq(y[i], y[ho])) {
                    ho = i;
                }
            }
            counter++
        }

        var V0 = math.subset(vertices, math.index(low, [0, size]));
        var y0 = f(this._toNormalVector(V0));
        return {vector: this._toNormalVector(V0), eval: y0};
    },

    /*
    * Function that get the local minimum of a function "f" with the gradient descent method
    *
    * Note: function based on the book Numerical Methods by John H. Mathews and Kurtis D. Fink,
    *       Spanish, Chap 8, pag 456.
    *
    * @param {Function} f - Function to get the local minimum, this function require one argument and should be a vector.
    * @param {Function} g - Function that represent the gradient of the function f
    * @param {Vector} p0 - Initial point approximation
    * @param {Number} maxIterations - Max iterations that we want
    * @param {Number|BigNumber} delta - Tolerance in process of search direction
    * @param {Number|BigNumber} epsilon - Tolerance for error y
    * @return {Object} result - return an object that contains a minimum vector (.vector) of the function and the evaluation of the function in this vector (.eval).
    *
    * */

    gradientDescent : function (f, g, p0, maxIterations, delta, epsilon) {
        var n = math.subset(math.size(p0), math.index(0));
        var maxj = 10;
        var big = 1e8;
        var h = 1;
        var P = math.zeros(maxj, n+1); // ojo con el tamaño
        var len = math.norm(p0);
        var y0 = f(p0);

        if(len > 1e4) h = math.divide(len, 1e4);

        var err = 1, counter = 0, condition = 0;
        math.subset(P, math.index(counter + 1, [0, n]), p0);
        math.subset(P, math.index(counter + 1, n), y0);

        while(counter < maxIterations && condition != 5 && (math.larger(h, delta) || math.larger(err, epsilon))) {
            var S = g(p0);

            var p1 = math.add(p0, math.multiply(h, S));
            var p2 = math.add(p0, math.multiply(math.multiply(h, 2), S));
            var y1 = f(p1);
            var y2 = f(p2);

            condition = 0;
            var j = 0;
            while(j < maxj && condition == 0) {
                len = math.norm(p0);

                if (math.smaller(y0, y1)) {
                    p2 = p1;
                    y2 = y1;
                    h = math.divide(h, 2);
                    p1 = math.add(p0, math.multiply(h, S));
                    y1 = f(p1);
                } else {
                    if(math.smaller(y2, y1)) {
                        p1 = p2;
                        y1 = y2;
                        h = math.multiply(2, h);
                        p2 = math.add(p0, math.multiply(math.multiply(2, h), S));
                        y2 = f(p2);
                    } else {
                        condition = -1;
                    }
                }
                j += 1;
                if (math.smaller(h, delta)) condition = 1;
                if (math.larger(math.abs(h), big) || math.larger(len, big)) condition = 5;
            }

            var pmin, ymin, hmin;

            if(condition == 5) {
                pmin = p1;
                ymin = y1;
            } else {
                var d = math.subtract(math.multiply(4, y1), math.subtract(math.multiply(2, y0), math.multiply(2, y2)));
                if(math.smaller(d, 0)) {
                    var mult = math.subtract(math.multiply(4, y1), math.subtract(math.multiply(3, y0), y2));
                    hmin = math.divide(math.multiply(h, mult), d);
                } else {
                    condition = 4;
                    hmin = math.divide(h, 3);
                }

                pmin = math.add(p0, math.multiply(hmin, S));
                ymin = f(pmin);

                var h0 = math.abs(hmin);
                var h1 = math.abs(math.subtract(hmin, h));
                var h2 = math.abs(math.subtract(hmin, math.multiply(2, h)));

                if(math.smaller(h0, h)) h = h0;
                if(math.smaller(h1, h)) h = h1;
                if(math.smaller(h2, h)) h = h2;
                if(math.equal(h, 0)) h = hmin;
                if(math.smaller(h, delta)) condition = 1;

                var e0 = math.abs(math.subtract(y0, ymin));
                var e1 = math.abs(math.subtract(y1, ymin));
                var e2 = math.abs(math.subtract(y2, ymin));

                if(math.unequal(e0, 0) && math.smaller(e0, err)) err = e0;
                if(math.unequal(e1, 0) && math.smaller(e1, err)) err = e1;
                if(math.unequal(e2, 0) && math.smaller(e2, err)) err = e2;
                if(math.equal(e0, 0) && math.equal(e1, 0) && math.equal(e2, 0)) err = 0;

                if(math.smaller(err, epsilon)) condition = 2;
                if(condition == 2 && math.smaller(h, delta)) condition = 3;
            }

            counter++;
            math.subset(P, math.index(counter + 1, [0, n]), pmin);
            math.subset(P, math.index(counter + 1, n), ymin);

            p0 = pmin;
            y0 = ymin;

            var ptemp = [];
            for (var i = 0; i < p0._data.length; ++i) ptemp.push(p0._data[i]);
            p0 = ptemp;
        }

        return {vector: p0, eval: y0};
    }
};
