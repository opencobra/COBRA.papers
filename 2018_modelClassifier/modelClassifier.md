
# Machine learning microbe classifier in Julia

This model classifier uses machine learning techniques to classifiy automatically microbe models.

## Import external and machine learning packages


```julia
Pkg.add("PyCall")
Pkg.add("DataFrames")
Pkg.add("CSV")
Pkg.add("ScikitLearn")
Pkg.add("PyPlot")
```


```julia
Pkg.update()
```


```julia
# import Julia packages
using DataFrames, CSV

# import ScikitLearn packages
using ScikitLearn, PyCall
using ScikitLearn.CrossValidation: train_test_split
import ScikitLearn: CrossValidation 
@sk_import model_selection: cross_val_score  
@sk_import preprocessing: (LabelEncoder, StandardScaler)
@sk_import metrics: accuracy_score 
@sk_import linear_model: LogisticRegression 
@sk_import ensemble: (RandomForestClassifier, AdaBoostClassifier)
@sk_import tree: DecisionTreeClassifier 
@sk_import datasets: (make_moons, make_circles, make_classification)
@sk_import neighbors: KNeighborsClassifier
@sk_import svm: SVC
@sk_import naive_bayes: GaussianNB
@sk_import discriminant_analysis: (LinearDiscriminantAnalysis, QuadraticDiscriminantAnalysis)
@sk_import metrics: confusion_matrix;
@pyimport matplotlib.colors as mplc
```

## Load the numerical characteristics


```julia
using Conda
Conda.add("matplotlib")
```


```julia
# load the data and the phylum labels
data = CSV.read("Supplementary_Table_S5.csv")
labels = data[:phylum]
labelsArr = data[:phylum]

head(data, 10)
```

## Curate and prepare the input data


```julia
# define the names of the categories
namesVect = names(data)

# create 'Others' category

# convert all labels to an array
labelsArr = convert(Array, labelsArr) 

# combine certain phyla to a special category
labelsArr = recode(labelsArr, "Thaumarchaeota"  => "Other",
                              "Crenarchaeota"   => "Other",
                              "Verrucomicrobia" => "Other",
                              "Spirochaetes"    => "Other",
                              "Cyanobacteria"   => "Other",
                              "Euryarchaeota"   => "Other",
                              "Synergistetes"   => "Other",
                              "Tenericutes"     => "Other",
                              "Fusobacteria"    => "Other",
                              "Planctomycetes"  => "Other");

#print out the various categories
levels(labelsArr)
#by(labelsArr, _:)
```


```julia
function selectAndCurateDate(data, numDataSet, topoDataSet, geneDataSet)
    # select certain categories
    # and normalize the data
    #=
    Note:

        Omitted categories are:
         1: name
         7: compSparsityRatio
         9: colDensityRel
        15: precisionEstimate
        16: estLevel
        42: phylum
        43: year
        44: gammain (topoDataSet)
        45: gammaout (topoDataSet)
        46: jaccard (geneDataSet)
    =#

    categories = []

    if numDataSet
        categories = [categories; [2:6...]; 8; [10:14...]; [17:41...]]
    end
    if topoDataSet
        categories = [categories; [44:45...]]
    end
    if geneDataSet
        categories = [categories; 46]
    end

    println(namesVect[categories])

    # select only the categories that should be used
    labelencoder = LabelEncoder() 

    for col in categories # 1:size(categories,1) #categories 
        data[col] = fit_transform!(labelencoder, data[col]) 
    end
    return data, categories
end
```

## Define the classification model


```julia
# define a classification model

function classification_model(model, data, predictors; testSize=0.33, randomState=100, cvK=10, printLevel=1) 
     y = convert(Array, labelsArr) 
     X = convert(Array, data[predictors]) 
     # split into train and test
     X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=testSize, random_state=randomState)
    
     # fit the model
     fit!(model, X_train, y_train)

     # make predictions on training set
     predictions = predict(model, X_train)
    
     # print accuracy 
     accuracy = accuracy_score(predictions, y_train)
    
     if printLevel > 0
         println("\nTraining score (accuracy): $(round(accuracy*100, 2))%")
     end

     # cvK-fold cross validation 
     cross_score = cross_val_score(model, X_train, y_train, cv=cvK, scoring="accuracy")
    
     # print cross_val_score 
     meanCS = mean(cross_score)
     if printLevel > 0
         println("Cross validation score: $(round(meanCS*100, 2))%") 
     end

     # return predictions 
     y_predictions_test = predict(model, X_test)
     testScore = accuracy_score(y_predictions_test, y_test)

     y_pred = predict(model, X);
    
     if printLevel > 0
         println("Test verification score: $(round(testScore*100, 2))%") 
     end
    
     return model, accuracy, meanCS, testScore, y_pred, y, predictions, y_train
end
```

## Determine the appropriate estimator

Source: http://scikit-learn.org/stable/tutorial/machine_learning_map/index.html


```julia
namesClf = ["Nearest Neighbors", "Linear SVM", "RBF SVM", "Decision Tree",
         "Random Forest", "AdaBoost", "Naive Bayes", "Linear Discriminant Analysis",
         "Logistic Regression"]

classifiers = [
    KNeighborsClassifier(5),
    SVC(kernel="linear", C=0.025),
    SVC(gamma=5, C=1),
    DecisionTreeClassifier(max_depth=5),
    RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1),
    AdaBoostClassifier(),
    GaussianNB(),
    LinearDiscriminantAnalysis(),
    LogisticRegression()];
```

**Conclusion on the most appropriate estimator**

According to the estimator decision graph, it is the **KNeighbors classifier** that is the most appropriate estimator to build a ML microbe classifier. The question is to determine the most appropriate number of neighbors and determine the optimal parameters for the classifier.

## Tuning the KNeighbors classifier


```julia
function chooseBestClassifier(categories, data)
    maxCS = 0
    maxAcc = 0
    maxTS = 0
    maxTestSize = 0.15
    maxY = 0
    maxYpred = 0
    maxPred = 0
    maxYtrain = 0
    # set the maximum testSize
    maxK = 2
    maxcvKi = 5
    for cvKi = maxcvKi:5:15
        for testSizej = maxTestSize:0.05:0.9
            for k = maxK:15
                model = KNeighborsClassifier(k)
                (model, accuracy, cross_score, test_score, y_pred, y, predictions, y_train) = classification_model(model, data, namesVect[categories], 
                                                                                testSize=testSizej, randomState=200, 
                                                                                cvK=cvKi, 
                                                                                printLevel=0)
                if test_score > maxTS
                    maxModel = model
                    maxPred = predictions
                    maxYtrain = y_train
                    maxY = y
                    maxYpred = y_pred
                    maxCS = cross_score
                    maxAcc = accuracy
                    maxTS = test_score
                    maxTestSize = testSizej
                    maxK = k
                    maxcvKi = cvKi
                end
            end
        end
    end
    println(" The best kNN estimator is: ")
    println("  > Accuracy: $(round(maxAcc*100, 2))")
    println("  > Test score: $(round(maxTS*100, 2))")
    println("  > Cross-validation score: $(round(maxCS*100, 2))")
    println("  > Number of neighbors N: $maxK")
    println("  > Test size P = $maxTestSize")
    println("  > Cross-validation k: $maxcvKi-fold")

    cnf_matrix = confusion_matrix(maxYtrain, maxPred)
    return cnf_matrix    
end
```


```julia
using PyPlot
"""
This function prints and plots the confusion matrix.
Normalization can be applied by setting `normalize=True`.
"""
function safeProp(num, den)
    if den == 0.
        prop = 0
    else
        prop = num / den
    end
    return prop
end
    
function plot_confusion_matrix(cm, classes)
    m, n = size(cm)
    extended_cm = zeros(m+1, n+1)
    extended_cm[1:m, 1:n] = cm
    imshow(cm, interpolation="nearest", cmap=ColorMap("Blues"))
    tick_params(axis="x", pad=30)
    colorbar(pad=0.15)
    tick_marks = 0:size(classes,1)-1
    xticks(tick_marks, classes, rotation=65)
    yticks(tick_marks, classes)
    text(size(cm, 2), size(cm, 1)+0.4 , "Accuracy",
                     horizontalalignment="center",
                     color="black", rotation=65)
        text(-1.1, size(cm, 1) , "Reliability",
                     horizontalalignment="center",
                     color="black")

    thresh = maximum(cm) / 2.
    for i in 1:size(cm, 1)
        for j in 1:size(cm, 2)
            if cm[i, j] > thresh 
                color="white"
            else
                color="black"
            end
            text(j-1, i-1, cm[i, j],
                     horizontalalignment="center",
                     color=color)

            prop = safeProp(cm[i, j], sum(cm))
            text(j-1, i-0.75, @sprintf("%.1f", prop*100)*"%",
                     horizontalalignment="center",
                     color=color)
        end
    end
    for i in 1:size(cm, 1)
        prop = safeProp(cm[i,i], sum(cm[i, :]))    
        text(size(cm, 2), i-1 , @sprintf("%.1f",prop*100)*"%",
                     horizontalalignment="center",
                     color="black")   
    end
    for j in 1:size(cm, 2)
            prop = safeProp(cm[j,j], sum(cm[:, j]))
        text(j-1, size(cm, 1) , @sprintf("%.1f",prop*100)*"%",
                     horizontalalignment="center",
                     color="black")
    end
        prop = safeProp(sum(diag(cm)), sum(cm))
    text(size(cm, 2), size(cm, 1) , @sprintf("%.1f",prop*100)*"%",
                     horizontalalignment="center",
                     color="black")
    tight_layout()
    ylabel("True label")
    xlabel("Predicted label")
end
```


```julia
function run(data, numDataSet, topoDataSet, geneDataSet)
    data_norm, selectedCategories = selectAndCurateDate(data, numDataSet, topoDataSet, geneDataSet)
    cnf_matrix = chooseBestClassifier(selectedCategories, data_norm)
                
    PyPlot.close_figs()
    plot_confusion_matrix(cnf_matrix, levels(labelsArr));
                
    filename="confusionMatrix"
    if numDataSet
        filename = filename * "_numDataSet"
    end
    if topoDataSet
        filename = filename * "_topoDataSet"
    end
    if geneDataSet
        filename = filename * "_geneDataSet"
    end
    filename = filename * ".pdf"
    savefig(filename)
end
```


```julia
for numDataSet in [true, false]
    for topoDataSet in [true, false]
         for geneDataSet in [true, false]
            if numDataSet || topoDataSet || geneDataSet
                run(data, numDataSet, topoDataSet, geneDataSet)
            end
        end
    end
end
```

## References: 

- https://www.analyticsvidhya.com/blog/2017/10/comprehensive-tutorial-learn-data-science-julia-from-scratch/
- https://kevinzakka.github.io/2016/07/13/k-nearest-neighbor/#more-on-k
- https://towardsdatascience.com/train-test-split-and-cross-validation-in-python-80b61beca4b6
