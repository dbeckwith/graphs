/**
 * Created by djbsn_000 on 5/2/2015.
 */

///<reference path="libs/underscore.js/underscore.d.ts"/>

interface PLU {
    P:Matrix;
    L:Matrix;
    U:Matrix;
}

class Matrix {
    public static byFunc(rows:number, cols:number, f:(i:number, j:number)=>number):Matrix {
        var values = [];
        for (var i = 0; i < rows; i++) {
            values.push([]);
            for (var j = 0; j < cols; j++) {
                values[i].push(f(i, j));
            }
        }
        return new Matrix(values);
    }

    public static zero(rows:number, cols:number):Matrix {
        return Matrix.byFunc(rows, cols, (i, j) => 0);
    }

    public static identity(size:number):Matrix {
        return Matrix.byFunc(size, size, (i, j) => +(i === j));
    }

    private _values:number[][];
    private valuesTrans:number[][];
    private _rows:number;
    private _cols:number;
    private square:boolean;
    private symmetric:boolean;

    /**
     * Creates a new Matrix.
     * @param values the values of the matrix, must be rectangular
     */
    public constructor(values:number[][]) {
        this._values = values.slice(0);
        this._rows = this._values.length;
        this._cols = 0;
        this._values.forEach((row:number[]) => {
            if (this._cols === 0) this._cols = row.length;
            else if (row.length !== this._cols) throw new Error('Matrix is not rectangular');
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

    /**
     * The two-dimensional array of values of this matrix.
     * @returns {number[][]}
     */
    public get values():number[][] {
        var vs:number[][] = [];
        for (var i = 0; i < this.rows; i++) {
            vs.push([]);
            for (var j = 0; j < this.cols; j++)
                vs[i].push(this._values[i][j]);
        }
        return vs;
    }

    /**
     * Creates a copy of this matrix.
     * @returns {Matrix}
     */
    public copy():Matrix {
        return Matrix.byFunc(this.rows, this.cols, (i, j) => this.get(i, j));
    }

    /**
     * Gets the value in the given row and column of this matrix.
     * @param row
     * @param col
     * @returns {number}
     */
    public get(row:number, col:number):number {
        if (row < 0 || row >= this.rows)
            throw new Error('Row number of out range.');
        if (col < 0 || col >= this.cols)
            throw new Error('Column number of out range.');
        return this._values[row][col];
    }

    /**
     * Gets the given row of this matrix.
     * @param row
     * @returns {number[]}
     */
    public getRow(row:number):number[] {
        if (row < 0 || row >= this.rows)
            throw new Error('Row number of out range.');
        return this._values[row];
    }

    /**
     * Gets the given column of this matrix.
     * @param col
     * @returns {number[]}
     */
    public getCol(col:number):number[] {
        if (col < 0 || col >= this.cols)
            throw new Error('Column number of out range.');
        return this.valuesTrans[col];
    }

    /**
     * The number of rows in this matrix.
     * @returns {number}
     */
    public get rows():number {
        return this._rows;
    }

    /**
     * The number of columns in this matrix.
     * @returns {number}
     */
    public get cols():number {
        return this._cols;
    }

    /**
     * Whether this matrix is square or not. A matrix is square if the number of rows is equal to the number of columns.
     * @returns {boolean}
     */
    public get isSquare():boolean {
        return this.square;
    }

    /**
     * Whether this matrix is symmetric or not. A matrix is symmetric if when reflected across its diagonal it is the
     * same.
     * @returns {boolean}
     */
    public get isSymmetric():boolean {
        return this.symmetric;
    }

    /**
     * Sets the value in the given row and column to the given number.
     * @param row
     * @param col
     * @param value
     */
    public set(row:number, col:number, value:number):void {
        if (row < 0 || row >= this.rows)
            throw new Error('Row number of out range.');
        if (col < 0 || col >= this.cols)
            throw new Error('Column number of out range.');
        this._values[row][col] = value;
    }

    /**
     * The transpose of this matrix, which is this matrix flipped along its diagonal.
     * @returns {Matrix}
     */
    public get transpose():Matrix {
        if (this.isSymmetric) return this.copy();
        return new Matrix(this.valuesTrans);
    }

    /**
     * Calculates the addition of this matrix and another matrix. The two matrices must have the same dimensions.
     * @param other another matrix to add to this one, must be the same dimensions
     * @returns {Matrix}
     */
    public add(other:Matrix):Matrix {
        if (this.rows !== other.rows && this.cols !== other.cols)
            throw new Error('Matrix dimensions do not agree');
        return Matrix.byFunc(this.rows, this.cols, (i, j) => this.get(i, j) + other.get(i, j));
    }

    /**
     * Calculates the subtraction of this matrix and another matrix. The two matrices must have the same dimensions.
     * @param other another matrix to subtract from this one, must be the same dimensions
     * @returns {Matrix}
     */
    public sub(other:Matrix):Matrix {
        if (this.rows !== other.rows && this.cols !== other.cols)
            throw new Error('Matrix dimensions do not agree');
        return Matrix.byFunc(this.rows, this.cols, (i, j) => this.get(i, j) - other.get(i, j));
    }

    /**
     * Calculates the multiplication of this matrix and a number or matrix. If the argument is a number, the resulting
     * matrix will be the same size as this matrix but will have every entry multiplied by the given number. If the
     * argument is another matrix, the result will be the matrix multiplication of this matrix and the argument. In this
     * case, the given matrix must have the same number of rows as this matrix has columns.
     * @param other a number to scale this matrix by, or a matrix to multiply this matrix by
     * @returns {Matrix}
     */
    public mult(other:number|Matrix):Matrix {
        if (typeof other === 'number') {
            if (other === 1)
                return this.copy();
            if (other === 0)
                return Matrix.zero(this.rows, this.cols);
            return Matrix.byFunc(this.rows, this.cols, (i, j) => other * this.get(i, j));
        }
        else {
            if (this.cols !== other.rows)
                throw new Error('Matrix dimensions do not agree');
            return Matrix.byFunc(this.rows, other.cols, (i, j) => {
                var sum = 0;
                for (var k = 0; k < this.cols; k++)
                    sum += this.get(i, k) * other.get(k, j);
                return sum;
            })
        }
    }

    /**
     * Calculates the division of this matrix and a number or matrix. If the argument is a number, the resulting
     * matrix will be the same size as this matrix but will have every entry divided by the given number. If the
     * argument is another matrix, the result will be the matrix multiplication of this matrix and the matrix inverse of
     * the argument. In this case, the given matrix must have the same number of rows as this matrix has columns.
     * @param other a number to scale this matrix by its inv, or a matrix to multiply this matrix by its inverse
     * @returns {Matrix}
     */
    public div(other:number|Matrix):Matrix {
        if (typeof other === 'number') {
            if (other === 1)
                return this.copy();
            return Matrix.byFunc(this.rows, this.cols, (i, j) => other / this.get(i, j));
        }
        else {
            if (this.cols !== other.rows)
                throw new Error('Matrix dimensions do not agree');
            return this.mult(other.inv);
        }
    }

    /**
     * Calculates the power of this matrix to the given integer. If the integer is positive, the result is this matrix
     * multiplied that many times by itself. If it is negative, the result is the inverse of that. If it is zero, the
     * identity is returned. This matrix must be square in order to have powers.
     * @param n an integer to raise this matrix to the power of
     * @returns {Matrix}
     */
    public pow(n:number):Matrix {
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
        var inv:boolean = n < 0;
        if (inv)
            n = -n;
        var m = this;
        for (var i = 0; i < n - 1; i++)
            m = m.mult(this);
        return inv ? m.inv : m;
    }

    private _minorMatrix:Matrix = null;

    /**
     * The matrix of minors of this matrix, where the i,j-th cell of the minor matrix is the i,j-th minor of
     * this matrix.
     * @returns {Matrix}
     */
    public get minorMatrix():Matrix {
        return this._minorMatrix ||
            (this._minorMatrix = Matrix.byFunc(this.rows, this.cols,
                (i, j) => Matrix.byFunc(this.rows - 1, this.cols - 1,
                    (i2, j2) => this.get(i2 < i ? i2 : i2 + 1, j2 < j ? j2 : j2 + 1)).det));
    }

    /**
     * Gets the i,j-th minor of this matrix. The i,j-th minor is the determinant of the matrix formed by
     * removing the ith row and jth column.
     * @param i
     * @param j
     * @returns {number}
     */
    public minor(i:number, j:number):number {
        return this.minorMatrix.get(i, j);
    }

    private _cofactorMatrix:Matrix = null;

    /**
     * The matrix of cofactors of this matrix, where the i,j-th cell of the cofactor matrix is the i,j-th cofactor of
     * this matrix.
     * @returns {Matrix}
     */
    public get cofactorMatrix():Matrix {
        return this._cofactorMatrix ||
            (this._cofactorMatrix = Matrix.byFunc(this.rows, this.cols,
                (i, j) => this.minor(i, j) * ((i + j) % 2 === 0 ? 1 : -1)));
    }

    /**
     * Gets the i,j-th cofactor of this matrix.
     * @param i
     * @param j
     * @returns {number}
     */
    public cofactor(i:number, j:number):number {
        return this.cofactorMatrix.get(i, j);
    }

    private PLU:PLU = null;
    private PLU_rowChanges:number = null;
    private PLU_uDiag:number[] = null;

    private PLUDecomp():PLU {
        if (!this.isSquare)
            throw new Error('Matrix is not square');
        if (this.PLU === null) {
            var n = this.rows;
            var A = this.copy();
            var rowPerm:number[] = _.range(n);
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
                if (max_r !== k) { // swap this row and that row in A and P
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
                P: Matrix.byFunc(n, n, (i, j) => +(rowPerm[j] === i)),
                L: Matrix.byFunc(n, n, (i, j) => i >= j ? (i === j ? 1 : A.get(i, j)) : 0),
                U: Matrix.byFunc(n, n, (i, j) => i <= j ? A.get(i, j) : 0)
            }
        }
        return this.PLU;
    }

    private _det:number = null;

    /**
     * The determinant of this matrix. This matrix must be square in order to have a det.
     * @returns {number}
     */
    public get det():number {
        if (!this.isSquare)
            throw new Error('Matrix is not square');
        if (this._det === null) {
            if (this.rows === 1)
                this._det = this.get(0, 0);
            if (this.rows === 2)
                this._det = this.get(0, 0) * this.get(1, 1) - this.get(0, 1) * this.get(1, 0);
            else if (this.rows === 3)
                this._det = this.get(0, 0) * this.get(1, 1) * this.get(2, 2) +
                    this.get(0, 1) * this.get(1, 2) * this.get(2, 0) +
                    this.get(0, 2) * this.get(1, 0) * this.get(2, 1) -
                    this.get(0, 2) * this.get(1, 1) * this.get(2, 0) -
                    this.get(0, 1) * this.get(1, 0) * this.get(2, 2) -
                    this.get(0, 0) * this.get(1, 2) * this.get(2, 1);
            else {
                if (this.rows <= 5) { // Laplace's formula, O(n!) is more efficient at n <= 5
                    console.log('using Laplace');
                    this._det = 0;
                    var i = 0;
                    for (var j = 0; j < this.rows; j++)
                        if (this.get(i, j) !== 0)
                            this._det += this.get(i, j) * this.cofactor(i, j);
                }
                else { // use PLU decomp O(n^3)
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
    }

    /**
     * Whether this matrix is invertible or not. A matrix is invertible if and only if it is square and its determinant
     * is non-zero.
     * @returns {boolean}
     */
    public get isInvertible():boolean {
        return this.isSquare && this.det !== 0;
    }

    private _adjugate:Matrix = null;

    /**
     * The adjugate matrix of this matrix, which is the transpose of the cofactor matrix.
     * @returns {Matrix}
     */
    public get adjugate():Matrix {
        return this._adjugate || (this._adjugate = this.cofactorMatrix.transpose);
    }

    private _inv:Matrix = null;

    /**
     * The inverse matrix of this matrix, which when multiplied by this matrix gives the identity. This matrix must be
     * square and have non-zero determinant in order to have an inverse.
     * @returns {Matrix}
     */
    public get inv():Matrix {
        if (!this.isInvertible)
            throw new Error('Matrix is non-invertible');
        return this._inv || (this._inv = this.adjugate.div(this.det));
    }

    private _trace:number = null;

    /**
     * The trace of this matrix. The trace of a matrix is the sum of its diagonal elements.
     * @returns {number}
     */
    public get trace():number {
        if (!this.square)
            throw new Error('Matrix is not square');
        if (this._trace === null) {
            this._trace = 0;
            for (var i = 0; i < this.rows; i++) {
                this._trace += this.get(i, i);
            }
        }
        return this._trace;
    }

    /**
     * Gets a plain-text string representation of this matrix in array form.
     * @returns {string}
     */
    public toString():string {
        return '[' + _.map(this.values, (row) => '[' + row.join(', ') + ']').join(', ') + ']';
    }

    /**
     * Gets a LaTeX-encoded string representation of this matrix.
     * @returns {string}
     */
    public toLatex():string {
        return '\\begin{bmatrix}' + _.map(this.values, (row) => row.join('&')).join('\\\\') + '\\end{bmatrix}';
    }

}
