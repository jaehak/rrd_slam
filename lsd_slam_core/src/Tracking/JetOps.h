#ifndef JETOPS_H
#define JETOPS_H

// JetOps from sample.h in libmv for accessing Jet<T,N> scalar value
// Added by Jae-Hak Kim, 2014-11-05
namespace ceres
{

// A jet traits class to make it easier to work with mixed auto / numeric diff.
template<typename T>
struct JetOps
{

    static bool IsScalar()
    {

        return true;

    }

    static T GetScalar(const T& t)
    {

        return t;

    }
    static void SetScalar(const T& scalar, T* t)
    {

        *t = scalar;

    }
    static void ScaleDerivative(double scale_by, T *value)
    {

        // For double, there is no derivative to scale.

    }

};

template<typename T, int N>
struct JetOps<Jet<T, N> >
{

    static bool IsScalar()
    {

        return false;

    }

    static T GetScalar(const Jet<T, N>& t)
    {

        return t.a;

    }

    static void SetScalar(const T& scalar, Jet<T, N>* t)
    {

        t->a = scalar;

    }

    static void ScaleDerivative(double scale_by, Jet<T, N> *value)
    {

        value->v *= scale_by;

    }

};

template<typename FunctionType, int kNumArgs, typename ArgumentType>
struct Chain
{

    static ArgumentType Rule(const FunctionType &f,
                             const FunctionType dfdx[kNumArgs],
                             const ArgumentType x[kNumArgs])
    {

        // In the default case of scalars, there's nothing to do since there
        // are no
        // derivatives to propagate.
        return f;

    }

};

// XXX Add documentation here!
template<typename FunctionType, int kNumArgs, typename T, int N>
struct Chain<FunctionType, kNumArgs, Jet<T, N> >
{

    static Jet<T, N> Rule(const FunctionType &f,
                          const FunctionType dfdx[kNumArgs],
                          const Jet<T, N> x[kNumArgs])
    {
        // x is itself a function of another variable ("z"); what this function
        // needs to return is "f", but with the derivative with respect to z
        // attached to the jet. So combine the derivative part of x's jets to
        // form
        // a Jacobian matrix between x and z (i.e. dx/dz).
        Eigen::Matrix<T, kNumArgs, N> dxdz;
        for (int i = 0; i < kNumArgs; ++i)
        {

            dxdz.row(i) = x[i].v.transpose();

        }

        // Map the input gradient dfdx into an Eigen row vector.
        Eigen::Map<const Eigen::Matrix<FunctionType, 1, kNumArgs> >
                vector_dfdx(dfdx, 1, kNumArgs);

        // Now apply the chain rule to obtain df/dz. Combine the derivative with
        // the scalar part to obtain f with full derivative information.
        Jet<T, N> jet_f;
        jet_f.a = f;
        jet_f.v = vector_dfdx.template cast<T>() * dxdz;  // Also known as dfdz.
        return jet_f;

    }

};

} // namespace ceres

#endif // JETOPS_H
