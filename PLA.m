[c, d]= perceptron_experiment(100, 10, 1000);
figure(1)
histogram(c);
xlabel('number of iterations');
ylabel('frequency');
title('Histogram of number of iterations for PLA in 1000 trails')
figure(2)
histogram(log(d));
xlabel('log of theoretical bound of interations minus actual number of iterations')
ylabel('frequency')
title('Histogram of log of theoretical bound of interactions minus actual number of interactions in 1000 trails')

function[w, iterations] = perceptron_learn(datain)
[nrow, ncol] = size(datain);
truelabel = datain(1:nrow, ncol);
truelabel = truelabel';
w = zeros(1, 11);
iterations = 0;
test = datain(1:nrow, 1:ncol-1)';
while isequal(sign(w*datain(1:nrow, 1:ncol-1)'), truelabel)==0
    index = find(sign(w*test)- truelabel ~= 0) ;
    w = w + datain(index(1), ncol )*datain(index(1), 1:ncol-1);
    iterations = iterations+1;
end
end

function[numiters, boundminusni] = perceptron_experiment(N, d, numsamples)
numiters = zeros(1, numsamples);
bounds = zeros(1, numsamples ) ;
for i = 1:numsamples
rng(i);
wstar = rand (1 ,d+1);
wstar(1)= 0;
x = 2*rand(d+1,N)-1 ;
x (1, 1:N)= ones(1,N) ;
correcttag = sign(wstar*x );
datain = [x', correcttag']; % nrow = N, ncol=d+2
[w, iter] = perceptron_learn(datain);
numiters(i)= iter ;
rho=min(correcttag.*(wstar*x)); %from problem 1 . 3
R2=max(sum(x.*x));
W2=sum(wstar.*wstar);
bounds (i)=R2*W2/rho ^ 2;
boundminusni = bounds-numiters;
end
end