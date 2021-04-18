import React, { Component } from "react";
import { Route, Switch, Redirect } from "react-router-dom";
import LoginForm from "./components/pages/login";
import MainPage from "./components/pages/mainpage";
import { ProtectedRoute } from "./components/common/protectedRoutes";
import Login1 from "./components/pages/login1";
import axios from "axios";

class App extends Component {
  render() {
    return (
      <div>
        <Switch>
          <Route path="/mainpage" component={MainPage} />
          <Route path="/login" component={Login1} />
          <Route path="/not-found" render={() => <h1>Not Found</h1>} />
          <Route
            path="/not-authorized"
            render={() => <h1>Not Authorized. Contact Admin</h1>}
          />
          <Redirect from="/" exact to="/mainpage" />
          <Redirect to="/not-found" />
        </Switch>
      </div>
    );
  }
}

export default App;
