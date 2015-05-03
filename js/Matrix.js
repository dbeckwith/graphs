/**
 * Created by djbsn_000 on 5/2/2015.
 */
///<reference path="libs/underscore.js/underscore.d.ts"/>
var Matrix = (function () {
    /**
     * Creates a new Matrix.
     * @param values the values of the matrix, must be rectangular
     */
    function Matrix(values) {
        var _this = this;
        this._minorMatrix = null;
        this._cofactorMatrix = null;
        this.PLU = null;
        this.PLU_rowChanges = null;
        this.PLU_uDiag = null;
        this._det = null;
        this._adjugate = null;
        this._inv = null;
        this._trace = null;
        this._values = values.slice(0);
        this._rows = this._values.length;
        this._cols = 0;
        this._values.forEach(function (row) {
            if (_this._cols === 0)
                _this._cols = row.length;
            else if (row.length !== _this._cols)
                throw new Error('Matrix is not rectangular');
        });
        this.square = this._rows === this._cols;
        this.valuesTrans = [];
        for (var j = 0; j < this._cols; j++) {
            this.valuesTrans.push([]);
            for (var i = 0; i < this._rows; i++) {
                this.valuesTrans[j].push(this._values[i][j]);
            }
        }
        this.symmetric = this.square;
        for (var i = 1; this.symmetric && i < this._rows; i++) {
            for (var j = 0; j < i; j++) {
                if (this.get(i, j) !== this.get(j, i)) {
                    this.symmetric = false;
                    break;
                }
            }
        }
    }
    Matrix.byFunc = function (rows, cols, f) {
        var values = [];
        for (var i = 0; i < rows; i++) {
            values.push([]);
            for (var j = 0; j < cols; j++) {
                values[i].push(f(i, j));
            }
        }
        return new Matrix(values);
    };
    Matrix.zero = function (rows, cols) {
        return Matrix.byFunc(rows, cols, function (i, j) { return 0; });
    };
    Matrix.identity = function (size) {
        return Matrix.byFunc(size, size, function (i, j) { return +(i === j); });
    };
    Object.defineProperty(Matrix.prototype, "values", {
        /**
         * The two-dimensional array of values of this matrix.
         * @returns {number[][]}
         */
        get: function () {
            var vs = [];
            for (var i = 0; i < this.rows; i++) {
                vs.push([]);
                for (var j = 0; j < this.cols; j++)
                    vs[i].push(this._values[i][j]);
            }
            return vs;
        },
        enumerable: true,
        configurable: true
    });
    /**
     * Creates a copy of this matrix.
     * @returns {Matrix}
     */
    Matrix.prototype.copy = function () {
        var _this = this;
        return Matrix.byFunc(this.rows, this.cols, function (i, j) { return _this.get(i, j); });
    };
    /**
     * Gets the value in the given row and column of this matrix.
     * @param row
     * @param col
     * @returns {number}
     */
    Matrix.prototype.get = function (row, col) {
        if (row < 0 || row >= this.rows)
            throw new Error('Row number of out range.');
        if (col < 0 || col >= this.cols)
            throw new Error('Column number of out range.');
        return this._values[row][col];
    };
    /**
     * Gets the given row of this matrix.
     * @param row
     * @returns {number[]}
     */
    Matrix.prototype.getRow = function (row) {
        if (row < 0 || row >= this.rows)
            throw new Error('Row number of out range.');
        return this._values[row];
    };
    /**
     * Gets the given column of this matrix.
     * @param col
     * @returns {number[]}
     */
    Matrix.prototype.getCol = function (col) {
        if (col < 0 || col >= this.cols)
            throw new Error('Column number of out range.');
        return this.valuesTrans[col];
    };
    Object.defineProperty(Matrix.prototype, "rows", {
        /**
         * The number of rows in this matrix.
         * @returns {number}
         */
        get: function () {
            return this._rows;
        },
        enumerable: true,
        configurable: true
    });
    Object.defineProperty(Matrix.prototype, "cols", {
        /**
         * The number of columns in this matrix.
         * @returns {number}
         */
        get: function () {
            return this._cols;
        },
        enumerable: true,
        configurable: true
    });
    Object.defineProperty(Matrix.prototype, "isSquare", {
        /**
         * Whether this matrix is square or not. A matrix is square if the number of rows is equal to the number of columns.
         * @returns {boolean}
         */
        get: function () {
            return this.square;
        },
        enumerable: true,
        configurable: true
    });
    Object.defineProperty(Matrix.prototype, "isSymmetric", {
        /**
         * Whether this matrix is symmetric or not. A matrix is symmetric if when reflected across its diagonal it is the
         * same.
         * @returns {boolean}
         */
        get: function () {
            return this.symmetric;
        },
        enumerable: true,
        configurable: true
    });
    /**
     * Sets the value in the given row and column to the given number.
     * @param row
     * @param col
     * @param value
     */
    Matrix.prototype.set = function (row, col, value) {
        if (row < 0 || row >= this.rows)
            throw new Error('Row number of out range.');
        if (col < 0 || col >= this.cols)
            throw new Error('Column number of out range.');
        this._values[row][col] = value;
    };
    Object.defineProperty(Matrix.prototype, "transpose", {
        /**
         * The transpose of this matrix, which is this matrix flipped along its diagonal.
         * @returns {Matrix}
         */
        get: function () {
            if (this.isSymmetric)
                return this.copy();
            return new Matrix(this.valuesTrans);
        },
        enumerable: true,
        configurable: true
    });
    /**
     * Calculates the addition of this matrix and another matrix. The two matrices must have the same dimensions.
     * @param other another matrix to add to this one, must be the same dimensions
     * @returns {Matrix}
     */
    Matrix.prototype.add = function (other) {
        var _this = this;
        if (this.rows !== other.rows && this.cols !== other.cols)
            throw new Error('Matrix dimensions do not agree');
        return Matrix.byFunc(this.rows, this.cols, function (i, j) { return _this.get(i, j) + other.get(i, j); });
    };
    /**
     * Calculates the subtraction of this matrix and another matrix. The two matrices must have the same dimensions.
     * @param other another matrix to subtract from this one, must be the same dimensions
     * @returns {Matrix}
     */
    Matrix.prototype.sub = function (other) {
        var _this = this;
        if (this.rows !== other.rows && this.cols !== other.cols)
            throw new Error('Matrix dimensions do not agree');
        return Matrix.byFunc(this.rows, this.cols, function (i, j) { return _this.get(i, j) - other.get(i, j); });
    };
    /**
     * Calculates the multiplication of this matrix and a number or matrix. If the argument is a number, the resulting
     * matrix will be the same size as this matrix but will have every entry multiplied by the given number. If the
     * argument is another matrix, the result will be the matrix multiplication of this matrix and the argument. In this
     * case, the given matrix must have the same number of rows as this matrix has columns.
     * @param other a number to scale this matrix by, or a matrix to multiply this matrix by
     * @returns {Matrix}
     */
    Matrix.prototype.mult = function (other) {
        var _this = this;
        if (typeof other === 'number') {
            if (other === 1)
                return this.copy();
            if (other === 0)
                return Matrix.zero(this.rows, this.cols);
            return Matrix.byFunc(this.rows, this.cols, function (i, j) { return other * _this.get(i, j); });
        }
        else {
            if (this.cols !== other.rows)
                throw new Error('Matrix dimensions do not agree');
            return Matrix.byFunc(this.rows, other.cols, function (i, j) {
                var sum = 0;
                for (var k = 0; k < _this.cols; k++)
                    sum += _this.get(i, k) * other.get(k, j);
                return sum;
            });
        }
    };
    /**
     * Calculates the division of this matrix and a number or matrix. If the argument is a number, the resulting
     * matrix will be the same size as this matrix but will have every entry divided by the given number. If the
     * argument is another matrix, the result will be the matrix multiplication of this matrix and the matrix inverse of
     * the argument. In this case, the given matrix must have the same number of rows as this matrix has columns.
     * @param other a number to scale this matrix by its inv, or a matrix to multiply this matrix by its inverse
     * @returns {Matrix}
     */
    Matrix.prototype.div = function (other) {
        var _this = this;
        if (typeof other === 'number') {
            if (other === 1)
                return this.copy();
            return Matrix.byFunc(this.rows, this.cols, function (i, j) { return other / _this.get(i, j); });
        }
        else {
            if (this.cols !== other.rows)
                throw new Error('Matrix dimensions do not agree');
            return this.mult(other.inv);
        }
    };
    /**
     * Calculates the power of this matrix to the given integer. If the integer is positive, the result is this matrix
     * multiplied that many times by itself. If it is negative, the result is the inverse of that. If it is zero, the
     * identity is returned. This matrix must be square in order to have powers.
     * @param n an integer to raise this matrix to the power of
     * @returns {Matrix}
     */
    Matrix.prototype.pow = function (n) {
        if (!this.isSquare)
            throw new Error('Matrix is not square');
        if (n % 1 !== 0)
            throw new Error('Only integer powers are supported');
        if (n === 0)
            return Matrix.identity(this.cols);
        if (n === 1)
            return this.copy();
        if (n === -1)
            return this.inv;
        var inv = n < 0;
        if (inv)
            n = -n;
        var m = this;
        for (var i = 0; i < n - 1; i++)
            m = m.mult(this);
        return inv ? m.inv : m;
    };
    Object.defineProperty(Matrix.prototype, "minorMatrix", {
        /**
         * The matrix of minors of this matrix, where the i,j-th cell of the minor matrix is the i,j-th minor of
         * this matrix.
         * @returns {Matrix}
         */
        get: function () {
            var _this = this;
            return this._minorMatrix || (this._minorMatrix = Matrix.byFunc(this.rows, this.cols, function (i, j) { return Matrix.byFunc(_this.rows - 1, _this.cols - 1, function (i2, j2) { return _this.get(i2 < i ? i2 : i2 + 1, j2 < j ? j2 : j2 + 1); }).det; }));
        },
        enumerable: true,
        configurable: true
    });
    /**
     * Gets the i,j-th minor of this matrix. The i,j-th minor is the determinant of the matrix formed by
     * removing the ith row and jth column.
     * @param i
     * @param j
     * @returns {number}
     */
    Matrix.prototype.minor = function (i, j) {
        return this.minorMatrix.get(i, j);
    };
    Object.defineProperty(Matrix.prototype, "cofactorMatrix", {
        /**
         * The matrix of cofactors of this matrix, where the i,j-th cell of the cofactor matrix is the i,j-th cofactor of
         * this matrix.
         * @returns {Matrix}
         */
        get: function () {
            var _this = this;
            return this._cofactorMatrix || (this._cofactorMatrix = Matrix.byFunc(this.rows, this.cols, function (i, j) { return _this.minor(i, j) * ((i + j) % 2 === 0 ? 1 : -1); }));
        },
        enumerable: true,
        configurable: true
    });
    /**
     * Gets the i,j-th cofactor of this matrix.
     * @param i
     * @param j
     * @returns {number}
     */
    Matrix.prototype.cofactor = function (i, j) {
        return this.cofactorMatrix.get(i, j);
    };
    Matrix.prototype.PLUDecomp = function () {
        if (!this.isSquare)
            throw new Error('Matrix is not square');
        if (this.PLU === null) {
            var n = this.rows;
            var A = this.copy();
            var rowPerm = _.range(n);
            this.PLU_uDiag = [];
            for (var i = 0; i < n; i++)
                this.PLU_uDiag.push(A.get(i, i));
            for (var k = 0; k < n - 1; k++) {
                var An = A.copy();
                // find max row in this column
                var max_r = k;
                for (var r = k + 1; r < n; r++) {
                    if (Math.abs(A.get(r, k)) > Math.abs(A.get(max_r, k)))
                        max_r = r;
                }
                if (max_r !== k) {
                    var temp = rowPerm[k];
                    rowPerm[k] = rowPerm[max_r];
                    rowPerm[max_r] = temp;
                    var temp2 = A._values[k];
                    A._values[k] = A._values[max_r];
                    A._values[max_r] = temp2;
                }
                for (var i = k + 1; i < n; i++) {
                    A.set(i, k, A.get(i, k) / A.get(k, k));
                    for (var j = k + 1; j < n; j++) {
                        An.set(i, j, A.get(i, j) - A.get(i, k) * A.get(k, j));
                        if (i === j)
                            this.PLU_uDiag[i] = An.get(i, j);
                    }
                }
                A = An.copy();
            }
            this.PLU_rowChanges = 0;
            for (var i = 0; i < n; i++)
                for (var j = i + 1; j < n; j++)
                    if (rowPerm[i] > rowPerm[j])
                        this.PLU_rowChanges++;
            this.PLU = {
                P: Matrix.byFunc(n, n, function (i, j) { return +(rowPerm[j] === i); }),
                L: Matrix.byFunc(n, n, function (i, j) { return i >= j ? (i === j ? 1 : A.get(i, j)) : 0; }),
                U: Matrix.byFunc(n, n, function (i, j) { return i <= j ? A.get(i, j) : 0; })
            };
        }
        return this.PLU;
    };
    Object.defineProperty(Matrix.prototype, "det", {
        /**
         * The determinant of this matrix. This matrix must be square in order to have a det.
         * @returns {number}
         */
        get: function () {
            if (!this.isSquare)
                throw new Error('Matrix is not square');
            if (this._det === null) {
                if (this.rows === 1)
                    this._det = this.get(0, 0);
                if (this.rows === 2)
                    this._det = this.get(0, 0) * this.get(1, 1) - this.get(0, 1) * this.get(1, 0);
                else if (this.rows === 3)
                    this._det = this.get(0, 0) * this.get(1, 1) * this.get(2, 2) + this.get(0, 1) * this.get(1, 2) * this.get(2, 0) + this.get(0, 2) * this.get(1, 0) * this.get(2, 1) - this.get(0, 2) * this.get(1, 1) * this.get(2, 0) - this.get(0, 1) * this.get(1, 0) * this.get(2, 2) - this.get(0, 0) * this.get(1, 2) * this.get(2, 1);
                else {
                    if (this.rows <= 5) {
                        console.log('using Laplace');
                        this._det = 0;
                        var i = 0;
                        for (var j = 0; j < this.rows; j++)
                            if (this.get(i, j) !== 0)
                                this._det += this.get(i, j) * this.cofactor(i, j);
                    }
                    else {
                        // (detA = det(P * L * U) = detP * detL * detU =
                        // (-1)^sgn(P) * prod l_ii * prod u_ii = (-1)^sgn(P) * prod u_ii)
                        console.log('using PLU');
                        this.PLUDecomp();
                        this._det = this.PLU_rowChanges % 2 === 0 ? 1 : -1;
                        for (var i = 0; this._det !== 0 && i < this.rows; i++) {
                            this._det *= this.PLU_uDiag[i];
                        }
                    }
                }
            }
            return this._det;
        },
        enumerable: true,
        configurable: true
    });
    Object.defineProperty(Matrix.prototype, "isInvertible", {
        /**
         * Whether this matrix is invertible or not. A matrix is invertible if and only if it is square and its determinant
         * is non-zero.
         * @returns {boolean}
         */
        get: function () {
            return this.isSquare && this.det !== 0;
        },
        enumerable: true,
        configurable: true
    });
    Object.defineProperty(Matrix.prototype, "adjugate", {
        /**
         * The adjugate matrix of this matrix, which is the transpose of the cofactor matrix.
         * @returns {Matrix}
         */
        get: function () {
            return this._adjugate || (this._adjugate = this.cofactorMatrix.transpose);
        },
        enumerable: true,
        configurable: true
    });
    Object.defineProperty(Matrix.prototype, "inv", {
        /**
         * The inverse matrix of this matrix, which when multiplied by this matrix gives the identity. This matrix must be
         * square and have non-zero determinant in order to have an inverse.
         * @returns {Matrix}
         */
        get: function () {
            if (!this.isInvertible)
                throw new Error('Matrix is non-invertible');
            return this._inv || (this._inv = this.adjugate.div(this.det));
        },
        enumerable: true,
        configurable: true
    });
    Object.defineProperty(Matrix.prototype, "trace", {
        /**
         * The trace of this matrix. The trace of a matrix is the sum of its diagonal elements.
         * @returns {number}
         */
        get: function () {
            if (!this.square)
                throw new Error('Matrix is not square');
            if (this._trace === null) {
                this._trace = 0;
                for (var i = 0; i < this.rows; i++) {
                    this._trace += this.get(i, i);
                }
            }
            return this._trace;
        },
        enumerable: true,
        configurable: true
    });
    /**
     * Gets a plain-text string representation of this matrix in array form.
     * @returns {string}
     */
    Matrix.prototype.toString = function () {
        return '[' + _.map(this.values, function (row) { return '[' + row.join(', ') + ']'; }).join(', ') + ']';
    };
    /**
     * Gets a LaTeX-encoded string representation of this matrix.
     * @returns {string}
     */
    Matrix.prototype.toLatex = function () {
        return '\\begin{bmatrix}' + _.map(this.values, function (row) { return row.join('&'); }).join('\\\\') + '\\end{bmatrix}';
    };
    return Matrix;
})();
//# sourceMappingURL=Matrix.js.map