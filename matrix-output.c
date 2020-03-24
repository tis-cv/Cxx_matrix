int main(void)
{
  int __retres;
  struct Matrix<2, 2> matrix_a;
  struct Matrix<2, 2> id;
  _Bool has_inverse;
  struct Matrix<2, 2> matrix_b;
  struct Matrix<2, 1> res;
  std::__1::ostream *tmp_3;
  struct std::__1::basic_ostream<char, std::__1::char_traits<char>> *tmp_2;
  _Bool tmp_4;
  __tis_globinit();
  Matrix<2, 2>::Ctor<double, double, double, double>(& matrix_a,2.,1.,4.,2.);
  identity<2>(& id);
  {
    struct Matrix<2, 2> __tis_arg;
    _Bool tmp;
    {
      {
        Matrix<2, 2>::Ctor(& __tis_arg,(struct Matrix<2, 2> const *)(& id));
        tmp = is_invertible<2>(& __tis_arg);
      }
      has_inverse = tmp;
    }
  }
  {
    char const *__tis_arg_0;
    struct std::__1::basic_ostream<char, std::__1::char_traits<char>> *tmp_1;
    struct std::__1::basic_ostream<char, std::__1::char_traits<char>> *tmp_0;
    {
      {
        if (has_inverse) __tis_arg_0 = "yes\n"; else __tis_arg_0 = "no\n";
      }
      tmp_0 = std::__1::operator<<<std::__1::char_traits<char>>(& std::__1::cout,
                                                                "identity is inversible: ");
    }
    tmp_1 = std::__1::operator<<<std::__1::char_traits<char>>(tmp_0,
                                                              __tis_arg_0);
  }
  {
    struct Matrix<2, 2> __tis_tmp_61;
    {
      {
        operator^(& __tis_tmp_61,(double)5,
                  (struct Matrix<2, 2> const *)(& id));
      }
      ;
      ;
    }
    operator+(& matrix_b,
              (struct Matrix_base<2, 2, Matrix> const *)(& matrix_a.__parent__Matrix_base<2, 2, Matrix>),
              (struct Matrix_base<2, 2, Matrix> const *)(& __tis_tmp_61.__parent__Matrix_base<2, 2, Matrix>));
  }
  {
    struct Matrix<2, 1> __tis_tmp_62;
    {
      { Matrix<2, 1>::Ctor<double, double>(& __tis_tmp_62,6.,10.); }
      ;
      ;
    }
    solve<2>(& res,(struct Matrix<2, 2> const *)(& matrix_b),
             (struct Matrix<2, 1> const *)(& __tis_tmp_62));
  }
  {
    ;
    tmp_2 = std::__1::operator<<<std::__1::char_traits<char>>(& std::__1::cout,
                                                              "RESULT IS:\n");
  }
  tmp_3 = operator<<(tmp_2,
                     (struct Matrix_base<2, 1, Matrix> const *)(& res.__parent__Matrix_base<2, 1, Matrix>));
  __retres = 0;
  goto return_label;
  tmp_4 = has_inverse;
  return_label: return __retres;
}