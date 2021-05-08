import React, { Component } from "react";
import { Route, Switch, Redirect } from "react-router-dom";
import LoginForm from "./components/pages/login";
import MainPage from "./components/pages/mainpage";
import { ProtectedRoute } from "./components/common/protectedRoutes";
import Login1 from "./components/pages/login1";
import axios from "axios";
import MainLayout from './components/common/main-layout'

class App extends Component {
  render() {
    return (
      <div className="App">
      <MainLayout></MainLayout>
    </div>
  )
  }
}

export default App;
