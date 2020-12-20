#Lab2 RS (V7)

setX = function(a, b, h){
  return(seq(from = a, to = b, by = h))
}

#Находим an, bn, cn, fn
calcAn = function(Ax, Bx, h){
  return(-(2 * Ax + h * Bx))
}

calcBn = function(Ax, Cx, h){
  return(4 * Ax + 2 * h * h * Cx)
}

calcCn = function(Ax, Bx, h){
  return(h * Bx - 2 * Ax)
}

calcFn = function(Fx, h){
  return(2 * h * h * Fx)
}

#Находим коэффциенты

#Нулевые
findAlpha1 = function(bn, cn){
  return(-cn[1]/bn[1])
}

findBeta1 = function(fn, bn){
  return(fn[1]/bn[1])
}


#Остальные
findK = function(an, bn, cn, fn, alpha, beta, N){
  alpha = vector("double", length = N)
  beta = vector("double", length = N)
  alpha[1] = findAlpha1(bn, cn)
  beta[1] = findBeta1(fn, bn)
  
  for(i in 2: N){
    alpha[i] = -cn[i]/(an[i] * alpha[i - 1] + bn[i])
    beta[i] = (fn[i] - an[i] * beta[i - 1])/(an[i] * alpha[i - 1] + bn[i])
  }
  return(data.frame(alpha, beta))
}

#Находим U до N включительно
calcU = function(an, bn, fn, alpha, beta, N){
  U = vector("double", length = N + 1);
  U[N + 1] = findUNN(an, bn, fn, alpha, beta, N)
  i = N
  while(i > 0){
    U[i] = alpha[i] * U[i + 1] + beta[i]
    i = i - 1
  }
  return(U)
}

#U[N + 1]:
findUNN = function(an, bn, fn, alpha, beta, N){
  return((fn[N + 1] - an[N + 1] * beta[N])/(an[N + 1] * alpha[N] + bn[N + 1]))
}

#РЕШЕНИЕ
solve = function(Ax, Bx, Cx, Fx, a, b, Ua, Ub, x){
  N = length(x) - 1
  h = (b - a)/(N)
  
  an = calcAn(Ax, Bx, h)
  bn = calcBn(Ax, Cx, h)
  cn = calcCn(Ax, Bx, h)
  fn = calcFn(Fx, h)
  
  #Задание начальных и конечных значений (для корректного поиска alpha и beta)
  an[N + 1] = 0
  bn[1] = 1
  bn[N + 1] = 1
  cn[1] = 0
  fn[1] = Ua
  fn[N + 1] = Ub
  
  #Нахождение коэффициентов
  tempFrame = findK(an, bn, cn, fn, alpha, beta, N)
  alpha = tempFrame[,1]
  beta = tempFrame[,2]
  
  #Нахождение значений U
  U = calcU(an, bn, fn, alpha, beta, N)
  
  #График
  plot(x, U, col = "white")
  plot(x, U)
  points(-0.5, U[1], col = "forestgreen")
  points(1.5, U[N + 1], col = "red")
  
  #Формирование датафрейма
  alphaN = c(alpha, NaN)
  betaN = c(beta, NaN)
  resFrame = data.frame(x, alphaN, betaN, U)
  return(resFrame)
}

#Решение: последовательность икс с шагом, решение и вывод датафрейма со всеми необходимыми данными
x = setX(-0.5, 1.5, 0.1)
solution = solve(exp(1 + sin(x)), cos(x), exp(-(2*x-1)*(2*x-1)/16), 1 - abs(cos(4*x)), -0.5, 1.5, 0.7, 1.1, x)
solution

#x1 = setX(-1, 1, 0.01)
#solution1 = solve(sqrt(2 + sin(2 * x)), 1 + sin(x), sqrt(3 + x), sin(x) * sin(x), -1, 1, 0.4, 0.8, x1)

