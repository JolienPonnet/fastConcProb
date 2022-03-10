// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

//////////////////////
// Count inversions //
//////////////////////
void merge(arma::vec& y, double nu, arma::uword left, arma::uword middle, arma::uword right, double& invCount, double& concCount, double& tiesObs, arma::vec& w, bool verb)
{
  // c is number of inversions
  arma::uword i, j, k, i2, j2, sTiesLeft, sTiesRight;
  // i2 points to the element in Left that is
  // at least nu larger than Right(j)
  // j2 points to the element in Right that is
  // at least nu larger than Left(i)
  //
  
  const arma::uword n = middle - left + 1; //length Left
  const arma::uword m = right - middle; // length Right
  arma::vec Left = y.subvec(left, left + n - 1);
  arma::vec Right = y.subvec(middle + 1, middle + m);
  
  arma::vec w_Left = w.subvec(left, left + n - 1);
  arma::vec w_Right = w.subvec(middle + 1, middle + m);
  
  i = 0, j = 0, k = left, i2 = 0, j2 = 0;
  sTiesLeft = 0, sTiesRight = 0;  // i is index in Left, j is index in Right, k is index in output of merge (y)
  
  //i2 = findFirstPos(Left.subvec(0, n-1), Right(0) + nu);
  //j2 = findFirstPos(Right.subvec(0, m-1), Left(0) + nu);
  
  double s_Left = arma::sum(w_Left);
  double s_Right = arma::sum(w_Right);
  
  // NOTE: i2 and j2 could also be initialized at 0
  // not sure which is faster, may depend on the data itself
  // it doesn't increase the complexity since there are
  // O(log(n)) merges and
  // findFirstPos (=binary search) is O(log(n))
  // so that adds O(log(n)^2) < O(nlog(n)) to the algorithm
  //
  while (i < n && j < m) // until end of Left or Right
  {
    // the loop first starts by making sure that
    // the counts i2 and j2 are correct.
    // note that neither i nor j are incremented there
    if(verb){
      Rcpp::Rcout << "i" << i <<  "\n";
      Rcpp::Rcout << "j" << j <<  "\n";
      Rcpp::Rcout << "k" << k <<  "\n";
      Rcpp::Rcout << "sTiesRight" << sTiesRight <<  "\n";
      Rcpp::Rcout << "sTiesRight" << sTiesRight <<  "\n";
      Rcpp::Rcout << "i2" << i2 <<  "\n";
      Rcpp::Rcout << "j2" << j2 <<  "\n";
      Rcpp::Rcout << "Left" << Left <<  "\n";
      Rcpp::Rcout << "Right" << Right <<  "\n";
    }
    
    if ((i2 < n && Left(i2) <= (Right(j) + nu))) { // no strict inequality, ties can also not be count as a discordant pair
      if(verb) Rcpp::Rcout << "update i2 " <<  "\n";
      if(Left(i2) - nu == Right(j)){
        sTiesLeft += w_Left(i2);
      }
      s_Left -= w_Left(i2);
      i2++;
    } else if ((j2 < m && Right(j2) <= (Left(i) + nu)) ) { // no strict inequality, ties can also not be count as a concordant pair
      if(verb) Rcpp::Rcout << "update j2 " <<   "\n";
      if(Right(j2) - nu == Left(i)){
        sTiesRight += w_Right(j2);
      }
      s_Right -= w_Right(j2);
      j2++;
    } else if (Left(i) <= Right(j)){ // Concordant
      if(verb) Rcpp::Rcout << "Left(i) <= Right(j) "  <<  "\n";
      //copy from the left and increase concordance count
      concCount+= w_Left(i)*s_Right;
      tiesObs+= w_Left(i)*sTiesRight;
  
      w(k) = w_Left(i);
      y(k++) = Left(i++);
      if(i < n && (Left(i) != Left(i-1))){
        sTiesRight = 0;
      }
    } else{ // Discordant
        if(verb) Rcpp::Rcout << "Left(i) > Right(j)"  <<  "\n";
      //copy from the right and increase inversion count
      invCount+= w_Right(j)*s_Left;
      if(nu > 0){// if nu == 0, tie is already counted in concordant part
        tiesObs+= w_Right(j)*sTiesLeft;
      }
      w(k) = w_Right(j);
      y(k++) = Right(j++);
      if(j < m && (Right(j) != Right(j-1))){
        sTiesLeft = 0;
      }
    }
    if(verb){
      Rcpp::Rcout << "concCount  "  << concCount<<  "\n";
      Rcpp::Rcout << "invCount "  << invCount<<  "\n";
      Rcpp::Rcout << "tiesObs "  << tiesObs<< "\n";
    }
  }
  
  while (i < n) // add remaining part of Left (already sorted)
  {
    if(verb) Rcpp::Rcout << "while i < n "  <<  "\n";
    w(k) = w_Left(i);
    y(k++) = Left(i++);
  }
  
  while (j < m) // add remaining part of Right (already sorted)
  {
    if(verb) Rcpp::Rcout << "while j < m "  <<  "\n";
    w(k) = w_Right(j);
    y(k++) = Right(j++);
  }
}

void mergeSort(arma::vec& y, double nu, arma::uword left, arma::uword right, double& invCount, double& concCount, double& tiesObs, arma::vec& w, bool verb)
{
  if (left < right)
  {
    arma::uword middle = left + (right - left) / 2;
    mergeSort(y, nu, left, middle, invCount, concCount, tiesObs, w, verb);
    mergeSort(y, nu, middle + 1, right, invCount, concCount, tiesObs,w, verb);
    merge(y, nu, left, middle, right, invCount, concCount, tiesObs, w, verb);
  }
}

// [[Rcpp::export]]
arma::vec countInversionsCpp(arma::vec y, arma::vec& w, double nu, bool verb)
{
  // calculate number of inversions in vector y
  // algorithm uses a modification of merge sort
  int n = y.size();
  
  arma::vec result = zeros(3);
  mergeSort(y, nu, 0, n - 1, result(0), result(1), result(2), w, verb); // index starts at 0
  return(result);
}






























// // [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>
// 
// using namespace arma;
// using namespace Rcpp;
// using namespace std;
// 
// //////////////////////
// // Count inversions //
// //////////////////////
// void merge(arma::vec& y, double nu, arma::uword left, arma::uword middle, arma::uword right, double& invCount, double& concCount, double& tiesObs, arma::vec& w, bool verb)
// {
//   // c is number of inversions
//   arma::uword i, j, k, i2, j2, sumTies, subInd, firstInd;
//   bool firstFound;
//   // i2 points to the element in Left that is
//   // at least nu larger than Right(j)
//   // j2 points to the element in Right that is
//   // at least nu larger than Left(i)
//   //
// 
//   const arma::uword n = middle - left + 1; //length Left
//   const arma::uword m = right - middle; // length Right
//   arma::vec Left = y.subvec(left, left + n - 1);
//   arma::vec Right = y.subvec(middle + 1, middle + m);
// 
//   arma::vec w_Left = w.subvec(left, left + n - 1);
//   arma::vec w_Right = w.subvec(middle + 1, middle + m);
// 
//   i = 0, j = 0, k = left, i2 = 0, j2 = 0;  // i is index in Left, j is index in Right, k is index in output of merge (y)
// 
//   //i2 = findFirstPos(Left.subvec(0, n-1), Right(0) + nu);
//   //j2 = findFirstPos(Right.subvec(0, m-1), Left(0) + nu);
// 
//   double s_Left = arma::sum(w_Left);
//   double s_Right = arma::sum(w_Right);
// 
//   // NOTE: i2 and j2 could also be initialized at 0
//   // not sure which is faster, may depend on the data itself
//   // it doesn't increase the complexity since there are
//   // O(log(n)) merges and
//   // findFirstPos (=binary search) is O(log(n))
//   // so that adds O(log(n)^2) < O(nlog(n)) to the algorithm
//   //
//   while (i < n && j < m) // tot je op einde bent van Left of Right
//   {
//     // the loop first starts by making sure that
//     // the counts i2 and j2 are correct.
//     // note that neither i nor j are incremented there
//     if(verb){
//       Rcpp::Rcout << "i" << i <<  "\n";
//       Rcpp::Rcout << "j" << j <<  "\n";
//       Rcpp::Rcout << "k" << k <<  "\n";
//       Rcpp::Rcout << "i2" << i2 <<  "\n";
//       Rcpp::Rcout << "j2" << j2 <<  "\n";
// 
//       Rcpp::Rcout << "Left" << Left <<  "\n";
//       Rcpp::Rcout << "Right" << Right <<  "\n";
//     }
// 
//     if ((i2 < n && Left(i2) <= (Right(j) + nu))) { // no longer strict inequality, ties can also not be count as a conconrdant or discordant pair
//       if(verb) Rcpp::Rcout << "update i2 " <<  "\n";
//       s_Left -= w_Left(i2);
//       i2++;
//     } else if ((j2 < m && Right(j2) <= (Left(i) + nu)) ) { // no longer strict inequality, ties can also not be count as a conconrdant or discordant pair
//       if(verb) Rcpp::Rcout << "update j2 " <<   "\n";
//       s_Right -= w_Right(j2);
//       j2++;
//     } else if(Left(i) == Right(j) + nu ){ // tie in observation
//         if(verb) Rcpp::Rcout << "Left(i) == Right(j) + nu "  <<  "\n";
//         // sumTies = sum(Right == Left(i)); // Since Right is sorted, while loop below is faster
//         subInd = i;
//         sumTies = 0;
//         while(subInd < n && Left(subInd) - nu == Right(j)){
//           subInd++;
//           sumTies++;
//         }
//         tiesObs+= w_Right(j)*sum(w_Left.subvec(i,(i+sumTies-1)));
//         //Right is smallest observation, so that one is selected
//         invCount+= w_Right(j)*s_Left;
//         w(k) = w_Right(j);
//         y(k++) = Right(j++);
//     } else if(Left(i) == Right(j) - nu){ // tie in observation
//         if(verb) Rcpp::Rcout << "Left(i) == Right(j) - nu "  <<  "\n";
//         // sumTies = sum(Right == Left(i)); // Since Left is sorted, while loop below is faster
//         subInd = j;
//         sumTies = 0;
//         while(subInd < m && Right(subInd) - nu == Left(i)){
//           subInd++;
//           sumTies++;
//         }
//         tiesObs+= w_Left(i)*sum(w_Right.subvec(j,(j+sumTies-1)));
//         //Left is smallest observation, so that one is selected
//         concCount+= w_Left(i)*s_Right;
//         w(k) = w_Left(i);
//         y(k++) = Left(i++);
//     } else if (Left(i) < Right(j)){ // Concordant
//       if(verb) Rcpp::Rcout << "Left(i) < Right(j) "  <<  "\n";
//       //copy from the left and increase concordance count
//       //concCount+= m - j2;
//       concCount+= w_Left(i)*s_Right;
//       w(k) = w_Left(i);
//       y(k++) = Left(i++);
//     } else if (Left(i) > Right(j)){ // Discordant
//       if(verb) Rcpp::Rcout << "Left(i) > Right(j)"  <<  "\n";
//       //copy from the right and increase inversion count
//       //invCount+= n - i2;
//       invCount+= w_Right(j)*s_Left;
//       w(k) = w_Right(j);
//       y(k++) = Right(j++);
//     }else{// Left(i) == Right(j) en nu > 0
//       if(verb) Rcpp::Rcout << "Left(i) == Right(j) "  <<  "\n";
//       // This comparison is no tie, but with definition of i2 en j2, don't forget any tie! e.g. merge 6, 11 and 6,8,10 with nu = 2. Left 6 is concordant with 10 and a tie with 8
//       subInd = j;
//       firstInd = j;
//       firstFound = FALSE;
//       sumTies = 0;
//       while(subInd < m && Right(subInd) - nu <= Left(i)){
//         if(Right(subInd) - nu == Left(i)){
//           sumTies++;
//           if(!firstFound){
//             firstFound = TRUE;
//             firstInd = subInd;
//           }
//         }
//         subInd++;
//       }
//       if(sumTies > 0){
//         tiesObs+= w_Left(i)*sum(w_Right.subvec(firstInd,(firstInd+sumTies-1)));
//       }
//       concCount+= w_Left(i)*s_Right;
//       w(k) = w_Left(i);
//       y(k++) = Left(i++);
//     }
//     if(verb){
//       Rcpp::Rcout << "concCount "  << concCount<<  "\n";
//       Rcpp::Rcout << "invCount "  << invCount<<  "\n";
//       Rcpp::Rcout << "tiesObs "  << tiesObs<< "\n";
//     }
//   }
// 
//   while (i < n) // add remaining part of Left (already sorted)
//   {
//     if(verb) Rcpp::Rcout << "while 1 "  <<  "\n";
//     w(k) = w_Left(i);
//     y(k++) = Left(i++);
//   }
// 
//   while (j < m) // add remaining part of Right (already sorted)
//   {
//     if(verb) Rcpp::Rcout << "while 2 "  <<  "\n";
//     w(k) = w_Right(j);
//     y(k++) = Right(j++);
//   }
// }
// 
// void mergeSort(arma::vec& y, double nu, arma::uword left, arma::uword right, double& invCount, double& concCount, double& tiesObs, arma::vec& w, bool verb)
// {
//   if (left < right)
//   {
//     arma::uword middle = left + (right - left) / 2;
//     mergeSort(y, nu, left, middle, invCount, concCount, tiesObs, w, verb);
//     mergeSort(y, nu, middle + 1, right, invCount, concCount, tiesObs,w, verb);
//     merge(y, nu, left, middle, right, invCount, concCount, tiesObs, w, verb);
//   }
// }
// 
// // [[Rcpp::export]]
// arma::vec countInversionsCpp(arma::vec y, arma::vec& w, double nu, bool verb)
// {
//   // calculate number of inversions in vector y
//   // algorithm uses a modification of merge sort
//   int n = y.size();
// 
//   arma::vec result = zeros(3);
//   mergeSort(y, nu, 0, n - 1, result(0), result(1), result(2), w, verb); // index starts at 0
//   return(result);
// }


// // [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>
// 
// using namespace arma;
// using namespace Rcpp;
// using namespace std;
// 
// //////////////////////
// // Count inversions //
// //////////////////////
// void merge(arma::vec& y, double nu, arma::uword left, arma::uword middle, arma::uword right, double& invCount, double& concCount, double& tiesObs, arma::vec& w, bool verb)
// {
//   // c is number of inversions
//   arma::uword i, j, k, i2, j2, nTiesLeft, nTiesRight, i2, j2, subIndLeft, subIndRight;
//   bool firstFoundLeft = FALSE;
//   bool firstFoundRight = FALSE;
//   // i2 points to the element in Left that is
//   // at least nu larger than Right(j)
//   // j2 points to the element in Right that is
//   // at least nu larger than Left(i)
//   //
//   
//   const arma::uword n = middle - left + 1; //length Left
//   const arma::uword m = right - middle; // length Right
//   arma::vec Left = y.subvec(left, left + n - 1);
//   arma::vec Right = y.subvec(middle + 1, middle + m);
//   
//   arma::vec w_Left = w.subvec(left, left + n - 1);
//   arma::vec w_Right = w.subvec(middle + 1, middle + m);
//   
//   i = 0, j = 0, k = left, i2 = 0, j2 = 0;
//   i2 = 0, j2 = 0, nTiesLeft = 0, nTiesRight = 0, subIndLeft = 0, subIndRight = 0;  // i is index in Left, j is index in Right, k is index in output of merge (y)
//   
//   //i2 = findFirstPos(Left.subvec(0, n-1), Right(0) + nu);
//   //j2 = findFirstPos(Right.subvec(0, m-1), Left(0) + nu);
//   
//   double s_Left = arma::sum(w_Left);
//   double s_Right = arma::sum(w_Right);
//   
//   // NOTE: i2 and j2 could also be initialized at 0
//   // not sure which is faster, may depend on the data itself
//   // it doesn't increase the complexity since there are
//   // O(log(n)) merges and
//   // findFirstPos (=binary search) is O(log(n))
//   // so that adds O(log(n)^2) < O(nlog(n)) to the algorithm
//   //
//   while (i < n && j < m) // tot je op einde bent van Left of Right
//   {
//     // the loop first starts by making sure that
//     // the counts i2 and j2 are correct.
//     // note that neither i nor j are incremented there
//     if(verb){
//       Rcpp::Rcout << "i" << i <<  "\n";
//       Rcpp::Rcout << "j" << j <<  "\n";
//       Rcpp::Rcout << "k" << k <<  "\n";
//       Rcpp::Rcout << "i2" << i2 <<  "\n";
//       Rcpp::Rcout << "j2" << j2 <<  "\n";
//       Rcpp::Rcout << "nTiesLeft" << nTiesLeft <<  "\n";
//       Rcpp::Rcout << "nTiesRight" << nTiesRight <<  "\n";
//       Rcpp::Rcout << "i2" << i2 <<  "\n";
//       Rcpp::Rcout << "j2" << j2 <<  "\n";
//       Rcpp::Rcout << "Left" << Left <<  "\n";
//       Rcpp::Rcout << "Right" << Right <<  "\n";
//     }
//     
//     if ((i2 < n && Left(i2) <= (Right(j) + nu))) { // no longer strict inequality, ties can also not be count as a conconrdant or discordant pair
//       if(verb) Rcpp::Rcout << "update i2 " <<  "\n";
//       if(Left(i2) - nu == Right(j)){
//         if(!firstFoundLeft){
//           firstFoundLeft = TRUE;
//           i2 = i2;
//         }
//         subIndLeft++;
//         nTiesLeft++;
//       }
//       s_Left -= w_Left(i2);
//       i2++;
//     } else if ((j2 < m && Right(j2) <= (Left(i) + nu)) ) { // no longer strict inequality, ties can also not be count as a conconrdant or discordant pair
//       if(verb) Rcpp::Rcout << "update j2 " <<   "\n";
//       if(Right(j2) - nu == Left(i)){
//         if(!firstFoundRight){
//           firstFoundRight = TRUE;
//           j2 = j2;
//         }
//         subIndRight++;
//         nTiesRight++;
//       }
//       s_Right -= w_Right(j2);
//       j2++;
//     } else if(Left(i) == Right(j) + nu ){ // tie in observation
//       if(verb) Rcpp::Rcout << "Left(i) == Right(j) + nu "  <<  "\n";
//       // sumTies = sum(Right == Left(i)); // Since Right is sorted, while loop below is faster
//       tiesObs+= w_Right(j)*sum(w_Left.subvec((i + i2),(i + i2 + nTiesLeft - 1)));
//       //Right is smallest observation, so that one is selected
//       invCount+= w_Right(j)*s_Left;
//       w(k) = w_Right(j);
//       y(k++) = Right(j++);
//       if(j < m && (Right(j) != Right(j-1))){
//         nTiesLeft = 0;
//         subIndLeft = 0;
//         firstFoundLeft = FALSE;
//         i2 = 0;
//       }
//     } else if(Left(i) == Right(j) - nu){ // tie in observation
//       if(verb) Rcpp::Rcout << "Left(i) == Right(j) - nu "  <<  "\n";
//       // sumTies = sum(Right == Left(i)); // Since Left is sorted, while loop below is faster
//       tiesObs+= w_Left(i)*sum(w_Right.subvec(j2,(j2 + nTiesRight - 1)));
//       //Left is smallest observation, so that one is selected
//       concCount+= w_Left(i)*s_Right;
//       w(k) = w_Left(i);
//       y(k++) = Left(i++);
//       if(i < n && (Left(i) != Left(i-1))){
//         nTiesRight = 0;
//         subIndRight = 0;
//         firstFoundRight = FALSE;
//         j2 = 0;
//     } else if (Left(i) < Right(j)){ // Concordant
//       if(verb) Rcpp::Rcout << "Left(i) < Right(j) "  <<  "\n";
//       //copy from the left and increase concordance count
//       //concCount+= m - j2;
//       concCount+= w_Left(i)*s_Right;
//       if(nTiesRight > 0){
//         tiesObs+= w_Left(i)*sum(w_Right.subvec(j2,(j2 + nTiesRight - 1)));
//       }
//       w(k) = w_Left(i);
//       y(k++) = Left(i++);
//       if(i < n && (Left(i) != Left(i-1))){
//         nTiesRight = 0;
//         subIndRight = 0;
//         firstFoundRight = FALSE;
//         j2 = 0;
//       }
//     } else if (Left(i) > Right(j)){ // Discordant
//       if(verb) Rcpp::Rcout << "Left(i) > Right(j)"  <<  "\n";
//       //copy from the right and increase inversion count
//       //invCount+= n - i2;
//       invCount+= w_Right(j)*s_Left;
//       if(nTiesLeft > 0 && nu > 0){// if nu == 0, tie is already counted in concordant part
//         tiesObs+= w_Right(j)*sum(w_Left.subvec(i2,(i2 + nTiesLeft - 1)));
//       }
//       w(k) = w_Right(j);
//       y(k++) = Right(j++);
//       if(j < m && (Right(j) != Right(j-1))){
//         nTiesLeft = 0;
//         subIndLeft = 0;
//         firstFoundLeft = FALSE;
//         i2 = 0;
//       }
//     } else{// Left(i) == Right(j) en nu > 0
//       if(verb) Rcpp::Rcout << "Left(i) == Right(j) "  <<  "\n";
//       // This comparison is no tie, but with definition of i2 en j2, don't forget any tie! e.g. merge 6, 11 and 6,8,10 with nu = 2. Left 6 is concordant with 10 and a tie with 8
//       if(nTiesRight > 0){
//         tiesObs+= w_Left(i)*sum(w_Right.subvec(j2,(j2 + nTiesRight - 1)));
//       }
//       concCount+= w_Left(i)*s_Right;
//       w(k) = w_Left(i);
//       y(k++) = Left(i++);
//       if(i < n && (Left(i) != Left(i-1))){
//         nTiesRight = 0;
//         subIndRight = 0;
//         firstFoundRight = FALSE;
//         j2 = 0;
//       }
//     }
//     if(verb){
//       Rcpp::Rcout << "concCount  "  << concCount<<  "\n";
//       Rcpp::Rcout << "invCount "  << invCount<<  "\n";
//       Rcpp::Rcout << "tiesObs "  << tiesObs<< "\n";
//     }
//   }
//   
//   while (i < n) // add remaining part of Left (already sorted)
//   {
//     if(verb) Rcpp::Rcout << "while 1 "  <<  "\n";
//     w(k) = w_Left(i);
//     y(k++) = Left(i++);
//   }
//   
//   while (j < m) // add remaining part of Right (already sorted)
//   {
//     if(verb) Rcpp::Rcout << "while 2 "  <<  "\n";
//     w(k) = w_Right(j);
//     y(k++) = Right(j++);
//   }
// }
// 
// void mergeSort(arma::vec& y, double nu, arma::uword left, arma::uword right, double& invCount, double& concCount, double& tiesObs, arma::vec& w, bool verb)
// {
//   if (left < right)
//   {
//     arma::uword middle = left + (right - left) / 2;
//     mergeSort(y, nu, left, middle, invCount, concCount, tiesObs, w, verb);
//     mergeSort(y, nu, middle + 1, right, invCount, concCount, tiesObs,w, verb);
//     merge(y, nu, left, middle, right, invCount, concCount, tiesObs, w, verb);
//   }
// }
// 
// // [[Rcpp::export]]
// arma::vec countInversionsCpp(arma::vec y, arma::vec& w, double nu, bool verb)
// {
//   // calculate number of inversions in vector y
//   // algorithm uses a modification of merge sort
//   int n = y.size();
//   
//   arma::vec result = zeros(3);
//   mergeSort(y, nu, 0, n - 1, result(0), result(1), result(2), w, verb); // index starts at 0
//   return(result);
// }
