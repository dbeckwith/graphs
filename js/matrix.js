
/* global _ */

function Matrix(values) {
  this.values = values;

  this.rows = values.length;
  this.cols = 0;
  values.forEach(function(row) {
    if (this.cols === 0) {
      this.cols = row.length;
    }
    else if (row.length !== this.cols) {
      throw new Error('Matrix is not rectangular');
    }
  }, this);
  this.square = this.rows === this.cols;
  this.valuesTrans = [];
  for (var j = 0; j < this.cols; j++) {
    this.valuesTrans.push([]);
    for (var i = 0; i < this.rows; i++) {
      this.valuesTrans[j].push(this.values[i][j]);
    }
  }

  this.copy = function() {
    return Matrix.byFunc(this.rows, this.cols, function(i, j) {
      return this.get(i, j);
    }, this);
  };

  this.get = function(row, col) {
    return this.values[row][col];
  };

  this.getRow = function(row) {
    return this.values[row];
  };

  this.getCol = function(col) {
    return this.valuesTrans[col];
  };

  this.symmetric = this.square;
  for (var i = 1; this.symmetric && i < this.rows; i++) {
    for (var j = 0; j < i; j++) {
      if (this.get(i, j) !== this.get(j, i)) {
        this.symmetric = false;
        break;
      }
    }
  }

  this.set = function(row, col, value) {
    this.values[row][col] = value;
  };

  this.transpose = function() {
    if (this.symmetric)
      return this;
    var a = this;
    return new Matrix(valuesTrans);
  };

  this.add = function(other) {
    if (this.rows !== other.rows && this.cols !== other.cols)
      throw new Error('Matrix dimensions do not agree');
    return Matrix.byFunc(this.rows, this.cols, function(i, j) {
      return this.get(i, j) + other.get(i, j);
    }, this);
  };

  this.sub = function(other) {
    if (this.rows !== other.rows && this.cols !== other.cols)
      throw new Error('Matrix dimensions do not agree');
    return Matrix.byFunc(this.rows, this.cols, function(i, j) {
      return this.get(i, j) - other.get(i, j);
    }, this);
  };

  this.mult = function(other) {
    if (other instanceof Number) {
      if (other === 1)
        return this;
      if (other === 0)
        return Matrix.zero(this.rows, this.cols);
      return Matrix.byFunc(this.rows, this.cols, function(i, j) {
        return other * this.get(i, j);
      }, this);
    }
    else if (other instanceof Matrix) {
      if (this.cols !== other.rows)
        throw new Error('Matrix dimensions do not agree');
      return Matrix.byFunc(this.rows, other.cols, function(i, j) {
        var sum = 0;
        for (var k = 0; k < this.cols; k++) {
          sum += this.get(i, k) * other.get(k, j);
        }
        return sum;
      }, this);
    }
  };

  this.div = function(other) {
    if (other instanceof Number) {
      if (other === 1)
        return this;
      return Matrix.byFunc(this.rows, this.cols, function(i, j) {
        return this.get(i, j) / other;
      }, this);
    }
    else if (other instanceof Matrix) {
      if (this.cols !== other.cols)
        throw new Error('Matrix dimensions do not agree');
      return this.mult(other.inverse());
    }
  };

  this.pow = function(n) {
    if (!this.square)
      throw new Error('Matrix is not square');
    if (n === 0)
      return Matrix.identity(this.cols);
    if (n === 1)
      return this;
    if (n === -1)
      return this.inverse();
    var inv = n < 0;
    if (inv)
      n = -n;
    var m = this;
    for (var i = 0; i < n - 1; i++) {
      m = m.mult(this);
    }
    return inv ? m.inverse() : m;
  };

  this.minor = function(i, j) {
    return Matrix.byFunc(this.rows - 1, this.cols - 1, function(i2, j2) {
      return this.get(i2 < i ? i2 : i2 + 1, j2 < j ? j2 : j2 + 1);
    }, this).determinant();
  };

  this.cofactor = function(i, j) {
    return this.minor(i, j) * ((i + j) % 2 === 0 ? 1 : -1);
  };

  var __PLU = null;
  var __PLU_rowChanges = null;
  var __PLU_uDiag = null;
  this.PLUDecomp = function() {
    if (!this.square)
      throw new Error('Matrix is not square');
    if (__PLU === null) {
      var n = this.rows;
      var A = this.copy();
      var rowPerm = _.range(n);
      __PLU_uDiag = [];
      for (var i = 0; i < n; i++) {
        __PLU_uDiag.push(A.get(i, i));
      }
      for (var k = 0; k < n - 1; k++) {
        var An = A.copy();
        // find max row in this column
        var max_r = k;
        for (var r = k + 1; r < n; r++) {
          if (Math.abs(A.get(r, k)) > Math.abs(A.get(max_r, k))) {
            max_r = r;
          }
        }
        if (max_r !== k) { // swap this row and that row in A and P
          var temp = rowPerm[k];
          rowPerm[k] = rowPerm[max_r];
          rowPerm[max_r] = temp;
          temp = A.values[k];
          A.values[k] = A.values[max_r];
          A.values[max_r] = A.values[k];
        }
        for (var i = k + 1; i < n; i++) {
          A.set(i, k, A.get(i, k) / A.get(k, k));
          for (var j = k + 1; j < n; j++) {
            An.set(i, j, A.get(i, j) - A.get(i, k) * A.get(k, j));
            if (i === j)
              __PLU_uDiag[i] = An.get(i, j);
          }
        }
        A = An.copy();
      }
      __PLU_rowChanges = 0;
      for (var i = 0; i < n; i++) {
        for (var j = i + 1; j < n; j++) {
          if (rowPerm[i] > rowPerm[j])
            __PLU_rowChanges++;
        }
      }
      __PLU = [
        Matrix.byFunc(n, n, function(i, j) {
          return +(rowPerm[j] === i);
        }),
        Matrix.byFunc(n, n, function(i, j) {
          return i >= j ? (i === j ? 1 : A.get(i, j)) : 0;
        }),
        Matrix.byFunc(n, n, function(i, j) {
          return i <= j ? A.get(i, j) : 0;
        })
      ];
    }
    return __PLU;
  };

  var __det = null;
  this.determinant = function() {
    if (!this.square)
      throw new Error('Matrix is not square');
    if (__det === null) {
      if (this.rows === 1)
        __det = this.get(0, 0);
      else if (this.rows === 2)
        __det = this.get(0, 0) * this.get(1, 1) - this.get(0, 1) * this.get(1, 0);
      else if (this.rows === 3)
        __det = this.get(0, 0) * this.get(1, 1) * this.get(2, 2) + this.get(0, 1) * this.get(1, 2) * this.get(2, 0) + this.get(0, 2) * this.get(1, 0) * this.get(2, 1) - this.get(0, 2) * this.get(1, 1) * this.get(2, 0) - this.get(0, 1) * this.get(1, 0) * this.get(2, 2) - this.get(0, 0) * this.get(1, 2) * this.get(2, 1);
      else {
        if (this.rows <= 5) { // Laplace's formula, O(n!) is more efficient at n <= 5
          console.log('using Laplace');
          __det = 0;
          var i = 0;
          for (var j = 0; j < this.rows; j++) {
            if (this.get(i, j) !== 0)
              __det += this.get(i, j) * this.cofactor(i, j);
          }
        }
        else { // use PLU decomp O(n^3) (detA = det(P * L * U) = detP * detL * detU = (-1)^sgn(P) * prod l_ii * prod u_ii = (-1)^sgn(P) * prod u_ii)
          console.log('using PLU');
          this.PLUDecomp();
          __det = __PLU_rowChanges % 2 === 0 ? 1 : -1;
          for (var i = 0; __det !== 0 && i < this.rows; i++) {
            __det *= __PLU_uDiag[i];
          }
        }
      }
    }
    return __det;
  };

  this.cofactorMatrix = function() {
    return Matrix.byFunc(this.rows, this.cols, function(i, j) {
      return this.cofactor(i, j);
    }, this);
  };

  this.adjugate = function() {
    return this.cofactorMatrix().transpose();
  };

  this.inverse = function() {
    return this.adjugate().div(this.determinant());
  };

  this.trace = function() {
    if (!this.square)
      throw new Error('Matrix is not square');
    var sum = 0;
    for (var i = 0; i < this.rows; i++) {
      sum += this.get(i, i);
    }
    return sum;
  };

  this.toString = function() {
    return '[' + _.map(this.values, function(row) {
      return '[' + row.join(', ') + ']';
    }).join(', ') + ']';
  };

  this.toLatex = function() {
    return '\\begin{bmatrix}' + _.map(this.values, function(row) {
      return row.join('&');
    }).join('\\\\') + '\\end{bmatrix}';
  };
}
Matrix.byFunc = function(rows, cols, f, context) {
  var values = [];
  for (var i = 0; i < rows; i++) {
    values.push([]);
    for (var j = 0; j < cols; j++) {
      values[i].push(f.apply(context, [i, j]));
    }
  }
  return new Matrix(values);
};
Matrix.zero = function(rows, cols) {
  return Matrix.byFunc(rows, cols, function(i, j) {
    return 0;
  });
};
Matrix.identity = function(size) {
  return Matrix.byFunc(size, size, function(i, j) {
    return +(i === j);
  });
};
