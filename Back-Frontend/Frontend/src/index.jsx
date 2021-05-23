import React from 'react';
import ReactDOM from 'react-dom';
import { BrowserRouter as Router, Route } from 'react-router-dom';
import { createLogger } from 'redux-logger';
import { createStore, applyMiddleware } from 'redux';
import { Provider } from 'react-redux';
import thunk from 'redux-thunk';
import reducers from '../src/store/reducers';
import { composeWithDevTools } from 'redux-devtools-extension';
import App from './App';
// import "./App.css";
import * as serviceWorker from './serviceWorker';

//bootstrap
// import "bootstrap/dist/css/bootstrap.css";
// import "bootstrap/dist/js/bootstrap";
// import $ from "jquery";
// import Popper from "popper.js";

//adminLTE
import 'jquery/src/jquery';
// import "admin-lte/bower_components/jquery/dist/jquery.min.js";
import 'admin-lte/bower_components/bootstrap/dist/js/bootstrap';
import 'admin-lte/bower_components/bootstrap/dist/css/bootstrap.css';
import 'admin-lte/bower_components/font-awesome/css/font-awesome.min.css';
import 'admin-lte/bower_components/Ionicons/css/ionicons.min.css';
import 'admin-lte/dist/css/AdminLTE.min.css';
import 'admin-lte/dist/css/skins/_all-skins.min.css';
import 'admin-lte/plugins/iCheck/icheck.min.js';
// import "admin-lte/bower_components/fastclick/lib/fastclick.js";
// import "admin-lte/dist/js/adminlte.min.js";

const middleware = [thunk];
if (process.env.NODE_ENV !== 'production') {
  middleware.push(createLogger());
}

const store = createStore(
  reducers,
  composeWithDevTools(applyMiddleware(...middleware)),
);

ReactDOM.render(
  <Provider store={store}>
    <Router>
      <App />
    </Router>
  </Provider>,
  document.getElementById('root'),
);

// If you want your app to work offline and load faster, you can change
// unregister() to register() below. Note this comes with some pitfalls.
// Learn more about service workers: https://bit.ly/CRA-PWA
serviceWorker.unregister();
