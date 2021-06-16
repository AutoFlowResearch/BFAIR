import React from 'react';
import { BrowserRouter as Router, Switch, Route, Link } from 'react-router-dom';
import DashboardPage from './components/templates/DashboardPage';
import InvestigationFormPage from './components/templates/InvestigationFormPage';

const Routes = () => (
  <Router>
    <Switch>
      <Route path='/doe'>
        <InvestigationFormPage />
      </Route>
      <Route path='/'>
        <DashboardPage />
      </Route>
      <Route path='/'></Route>
    </Switch>
  </Router>
);

export default Routes;
