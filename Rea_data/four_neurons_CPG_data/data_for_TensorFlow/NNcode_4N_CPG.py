import tensorflow as tf

#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

from sklearn import model_selection
import numpy as np
import csv

N = 200000 #how many points to use in training
numOfTestSampl = 28000

inputs_indicies = [0,1,2,3,4,5,6,7,8,9,10,11,12,13]
outputs_indicies = [14]

# Training Data:
dataMatrix = csv.reader(open('dataForTensorFlow_4N_trainValid_05_07_13_48.csv'), delimiter=",")
dataMatrix = list(dataMatrix)
dataMatrix = np.array(dataMatrix).astype("float")

sampl = dataMatrix[inputs_indicies,0:N] #0:27491
targ = dataMatrix[outputs_indicies,0:N]
sampl = np.transpose(sampl)
targ = np.transpose(targ)

# Test Data: (the data to use to calculate the MSE in the end)
dataMatrix_test = csv.reader(open('dataForTensorFlow_4N_test05_07_13_48.csv'), delimiter=",")
dataMatrix_test = list(dataMatrix_test)
dataMatrix_test = np.array(dataMatrix_test).astype("float")

sampl_test = np.transpose(dataMatrix_test[inputs_indicies,0:numOfTestSampl])
targ_test = np.transpose(dataMatrix_test[outputs_indicies,0:numOfTestSampl])


numIn = 14 #???

x, y = sampl, targ
test2trainRatio = 0.2 # in tensorflow this will be use for validation
X_train, X_test, Y_train, Y_test = model_selection.train_test_split(x, y, test_size=test2trainRatio)

total_len = X_train.shape[0]

# Parameters
training_epochs = 100
batch_size = 128
display_step = 5

# initialize error vectors:
train_err_vec = []
test_err_vec = []
epochNum = np.arange(training_epochs)+1

# tf Graph input
x = tf.placeholder(tf.float32, [None, numIn])
#g_ = tf.placeholder(tf.float32, [1,None])
y_ = tf.placeholder(tf.float32, [None, 1])

if False: # NN with 3 layers
	hidden = [10, 10, 10]
	out = 1

	W1 = tf.Variable(tf.random_uniform([numIn,hidden[0]],-1/numIn,1/numIn), name = 'net1/w1')
	b1 = tf.Variable(tf.random_uniform([hidden[0]],-1/numIn,1/numIn), name='net1/b1')
	z1 = tf.add(tf.matmul(x,W1), b1, name='net1/z1')
	h1 = tf.nn.relu(z1, name='net1/h1')

	W2 = tf.Variable(tf.random_uniform([hidden[0],hidden[1]],-1/hidden[0],1/hidden[0]),name = 'net1/w2')
	b2 = tf.Variable(tf.random_uniform([hidden[1]],-1/hidden[0],1/hidden[0]), name='net1/b2')
	z2 = tf.add(tf.matmul(h1,W2), b2,name='net1/z2')
	h2 = tf.nn.relu(z2, name='net1/h2')
	
	W3 = tf.Variable(tf.random_uniform([hidden[1],hidden[2]],-1/hidden[1],1/hidden[1]),name = 'net1/w3')
	b3 = tf.Variable(tf.random_uniform([hidden[2]],-1/hidden[1],1/hidden[1]), name='net1/b3')
	z3 = tf.add(tf.matmul(h2,W3), b3,name='net1/z3')
	h3 = tf.nn.relu(z3, name='net1/h2')
	
	W_out = tf.Variable(tf.random_uniform([hidden[2],out],-1/hidden[2],1/hidden[2]), name = 'net1/w_out')
	b_out = tf.Variable(tf.random_uniform([out],-1/hidden[2],1/hidden[2]), name='net1/b_out')
	z_out = tf.add(tf.matmul(h3,W_out), b_out, name='net1/z_out')
	
if False: # NN with 2 layers
	hidden = [10, 10]
	out = 1

	W1 = tf.Variable(tf.random_uniform([numIn,hidden[0]],-1/numIn,1/numIn), name = 'net1/w1')
	b1 = tf.Variable(tf.random_uniform([hidden[0]],-1/numIn,1/numIn), name='net1/b1')
	z1 = tf.add(tf.matmul(x,W1), b1, name='net1/z1')
	h1 = tf.nn.relu(z1, name='net1/h1')

	W2 = tf.Variable(tf.random_uniform([hidden[0],hidden[1]],-1/hidden[0],1/hidden[0]),name = 'net1/w2')
	b2 = tf.Variable(tf.random_uniform([hidden[1]],-1/hidden[0],1/hidden[0]), name='net1/b2')
	z2 = tf.add(tf.matmul(h1,W2), b2,name='net1/z2')
	h2 = tf.nn.relu(z2, name='net1/h2')
	
	W_out = tf.Variable(tf.random_uniform([hidden[1],out],-1/hidden[1],1/hidden[1]), name = 'net1/w_out')
	b_out = tf.Variable(tf.random_uniform([out],-1/hidden[1],1/hidden[1]), name='net1/b_out')
	z_out = tf.add(tf.matmul(h2,W_out), b_out, name='net1/z_out')

if True: # NN with 1 hidden layers
	hidden = 500
	out = 1

	W1 = tf.Variable(tf.random_uniform([numIn,hidden],-1/numIn,1/numIn), name = 'net1/w1')
	b1 = tf.Variable(tf.random_uniform([hidden],-1/numIn,1/numIn), name='net1/b1')
	z1 = tf.add(tf.matmul(x,W1), b1, name='net1/z1')
	h1 = tf.nn.relu(z1, name='net1/h1')
		
	W_out = tf.Variable(tf.random_uniform([hidden,out],-1/hidden,1/hidden), name = 'net1/w_out')
	b_out = tf.Variable(tf.random_uniform([out],-1/hidden,1/hidden), name='net1/b_out')
	z_out = tf.add(tf.matmul(h1,W_out), b_out, name='net1/z_out')

y = z_out

# ## regular NN:
Loss1 = tf.reduce_mean(tf.square(y-y_), name = 'Loss1_MSE')
# ## Loss function for MoE:
# Loss2 = -tf.log(tf.reduce_sum(1*tf.exp(-1/2*(tf.square(y-y_)))), name='Likelihood_loss') # need to add g instead of 1

# ## minimization Algorhitms:
optimizer1 = tf.train.GradientDescentOptimizer(learning_rate=0.01).minimize(Loss1)
#optimizer1 = tf.train.AdamOptimizer(learning_rate=0.001,beta1=0.9, beta2=0.999, epsilon=1e-07).minimize(Loss1)
# optimizer2 = tf.train.GradientDescentOptimizer(learning_rate=learning_rate).minimize(Loss2)

saver = tf.train.Saver()
# Launch the graph
with tf.Session() as sess:
	sess.run(tf.global_variables_initializer())
	
	# # save log file
	# tflogger = tf.summary.FileWriter("C:/Users/user/Dropbox/Importqant stuff for thesis/Work for Jonathans paper/python work/logs", sess.graph)

    # Training cycle
	for epoch in range(training_epochs):
		train_err = 0.
		total_batch = int(total_len/batch_size)
		# np.random.permutation(X_train)
		# np.random.permutation(Y_train)
        # # Loop over all batches
		for i in range(total_batch-1):
			batch_x = X_train[i*batch_size:(i+1)*batch_size]
			batch_y = Y_train[i*batch_size:(i+1)*batch_size]
			#g = np.ones(1,batch_x.shape[0])
			# Run optimization op (backprop) and cost op (to get loss value)
			Error = sess.run(Loss1, {x:batch_x, y_:batch_y })
			sess.run(optimizer1, {x:batch_x, y_:batch_y})
				
			# Compute average loss
			train_err += Error
				
		train_err /= i
		test_err = sess.run(Loss1, {x:X_test, y_:Y_test })
		
		train_err_vec = np.append(train_err_vec,train_err)
		test_err_vec = np.append(test_err_vec,test_err)		


		# Display logs per epoch step
		if epoch % display_step == 0:
			print("Epoch:", (epoch+1), "train_err=", \
				"{:.9f}".format(train_err))
			print("[*]----------------------------")
			print("test err=", "{:.9f}".format(test_err))
			print("[*]============================")
	print ("Optimization Finished!")
	
	# # save model file:
	# saver.save(sess, "C:/Users/user/Dropbox/Importqant stuff for thesis/Work for Jonathans paper/python work/logs/model.ckpt")

	outTest = sess.run(y,{x:sampl_test})
	targTest = targ_test
	
	# print('the len of outTest: '+str(outTest.shape))
	# print('the len of targTest: '+str(targTest.shape))
	
	plt.scatter(targTest,outTest,label='Regression graph')
	plt.ylabel('output')
	plt.xlabel('targets')
	plt.legend()
	plt.show(block=True)
    
	plt.plot(epochNum,train_err_vec,label='train err')
	plt.plot(epochNum,test_err_vec,label='test err')
	plt.ylabel('error MSE')
	plt.xlabel('epoch num')
	plt.legend()
	plt.show(block=True)

	
    # Test model
    #Err = sess.run(Loss1, {x:X_test, y_:Y_test}))
    #correct_prediction = tf.equal(tf.argmax(pred, 1), tf.argmax(y, 1))
    # Calculate accuracy
    #accuracy = tf.reduce_mean(tf.cast(correct_prediction, "float"))
    #print ("Accuracy:", Err) #accuracy.eval({x: X_test, y: Y_test}))